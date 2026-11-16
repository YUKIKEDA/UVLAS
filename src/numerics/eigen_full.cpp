#include "eigen_full_detail.hpp"
#include "eigen_full_phase_private.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <utility>

namespace uvlas::numerics::eigen {
namespace {
using Complex = detail::FullComplex;
using Matrix = detail::FullMatrix;
using Data = detail::FullData;
using core::FailureCode;
constexpr auto operation = core::Operation::ExtractEigen;
struct Abort {
    core::Failure failure;
    AttemptDiagnostics diagnostics{};
};
[[noreturn]] void reject(FailureCode code) {
    throw Abort{core::Failure{code, operation}};
}
template <class T>
T need(core::Result<T> r) {
    if (!r)
        throw Abort{r.error()};
    if constexpr (!std::is_void_v<T>)
        return std::move(*r);
}
template <class T>
T need(Result<T> r) {
    if (!r)
        throw Abort{r.error().cause, r.error().diagnostics};
    if constexpr (!std::is_void_v<T>)
        return std::move(*r);
}
template <class T, class F>
Result<T> boundary(Context c, StateDiagnosticContext s, AttemptDiagnostics& d, F&& f) noexcept {
    try {
        return f();
    } catch (const Abort& e) {
        return std::unexpected(complete_failure(EigenFailure{e.failure, e.diagnostics}, c, s, d));
    } catch (const std::bad_alloc&) {
        return std::unexpected(
            complete_failure(core::Failure{FailureCode::ResourceExceeded, operation}, c, s, d));
    } catch (...) {
        return std::unexpected(
            complete_failure(core::Failure{FailureCode::FactorizationFailed, operation}, c, s, d));
    }
}
struct ReportSync {
    AttemptDiagnostics& outer;
    const Data& data;
    ~ReportSync() {
        outer = data.diagnostics;
    }
};
using FullClock = std::chrono::steady_clock;
void add_timing(Data& d, core::StageTiming& timing, FullClock::time_point start) noexcept {
    if (d.clock_reads == std::numeric_limits<std::uint64_t>::max()) {
        d.timing_overflow = true;
        return;
    }
    ++d.clock_reads;
    const auto duration =
        std::chrono::duration_cast<std::chrono::nanoseconds>(FullClock::now() - start).count();
    if (duration < 0 || timing.calls == std::numeric_limits<std::uint64_t>::max()) {
        d.timing_overflow = true;
        return;
    }
    const auto elapsed = static_cast<std::uint64_t>(duration);
    if (elapsed > std::numeric_limits<std::uint64_t>::max() - timing.elapsed_nanoseconds) {
        d.timing_overflow = true;
        return;
    }
    ++timing.calls;
    timing.elapsed_nanoseconds += elapsed;
}
class FullPreparationClock {
public:
    explicit FullPreparationClock(FullTimingMode mode) {
        if (mode == FullTimingMode::Enabled)
            start_ = FullClock::now();
    }
    void finish(Data& d) const {
        if (start_ && d.timings) {
            d.clock_reads = 1;
            add_timing(d, d.timings->preparation, *start_);
        }
    }

private:
    std::optional<FullClock::time_point> start_;
};
class FullStageClock {
public:
    FullStageClock(Data& d, core::StageTiming FullTimings::* stage) : data_(d), stage_(stage) {
        if (d.timings) {
            if (d.clock_reads == std::numeric_limits<std::uint64_t>::max()) {
                d.timing_overflow = true;
                return;
            }
            ++d.clock_reads;
            start_ = FullClock::now();
        }
    }
    ~FullStageClock() {
        if (start_)
            add_timing(data_, (*data_.timings).*stage_, *start_);
    }
    FullStageClock(const FullStageClock&) = delete;
    FullStageClock& operator=(const FullStageClock&) = delete;

private:
    Data& data_;
    core::StageTiming FullTimings::* stage_;
    std::optional<FullClock::time_point> start_;
};
void check_timing(const Data& d) {
    if (d.timing_overflow)
        reject(FailureCode::ResourceExceeded);
}
bool finite(Complex z) {
    return std::isfinite(z.real()) && std::isfinite(z.imag());
}
void finite_value(double x) {
    if (!std::isfinite(x))
        reject(FailureCode::Nonfinite);
}
void requirement(MetricRequirement r) {
    finite_value(r.budget);
    if (r.budget < 0)
        reject(FailureCode::InvalidInput);
}
void backward_requirement(const BackwardRequirement& r) {
    if (r.valueless_by_exception())
        reject(FailureCode::InvalidInput);
    if (const auto* q = std::get_if<MetricRequirement>(&r))
        requirement(*q);
}
void association(FullAssociation a) {
    finite_value(a.absolute);
    finite_value(a.relative);
    if (a.absolute < 0 || a.relative < 0)
        reject(FailureCode::InvalidInput);
}
std::pair<core::NumericalGeneration, core::NumericalGeneration>
generations(const FullOutputRequest& r) {
    if (const auto* keep = std::get_if<linear::RetainedPrecision>(&r.returned))
        return {keep->generation, keep->generation};
    if (const auto* p = std::get_if<linear::Binary64Projection>(&r.returned)) {
        if (p->retained_generation == p->target_generation)
            reject(FailureCode::InvalidInput);
        return {p->retained_generation, p->target_generation};
    }
    reject(FailureCode::PrecisionUnsupported);
}
void output_preflight(core::Index n, Context c, const FullOutputRequest& r, core::Index count) {
    finite_value(c.time);
    if (r.timing != FullTimingMode::Disabled && r.timing != FullTimingMode::Enabled)
        reject(FailureCode::InvalidInput);
    if (c.cache.representation != core::Precision::Binary64)
        reject(FailureCode::PrecisionUnsupported);
    if (n < 0 || r.limits.maximum_order < 0 || r.limits.maximum_results < 0)
        reject(FailureCode::InvalidInput);
    if (n > r.limits.maximum_order || count > r.limits.maximum_results)
        reject(FailureCode::ResourceExceeded);
    const std::array dims{n, core::Index{30}};
    need(core::checked_extent_product(dims, operation));
    for (const auto& q : r.consumers)
        requirement(q);
    if (count == 0 && !r.consumers.empty())
        reject(FailureCode::PrecisionUnresolved);
    if (!r.consumers.empty() && !r.evaluator)
        reject(FailureCode::InvalidInput);
    if (count > 0 && !r.consumers.empty()) {
        const auto supported =
            require_precision(core::Precision::Binary64, r.evaluator->supported_precision());
        if (!supported)
            throw Abort{core::Failure{supported.error().code, operation}};
    }
    if (r.require_physical_coverage)
        reject(FailureCode::EigenCoverageUnresolved);
    (void)generations(r);
}
void preflight(core::Index n, Context c, const FullRequest& r, core::Index count) {
    output_preflight(n, c, r, count);
    requirement(r.basis.residual);
    requirement(r.basis.orthogonality);
    backward_requirement(r.basis.backward);
    if (r.basis.residual.budget > 1e-10)
        reject(FailureCode::InvalidInput);
}
core::Index square(core::Index n) {
    const std::array dims{n, n};
    return need(core::checked_extent_product(dims, operation));
}
void scale_inputs(core::Index n, FullScales s) {
    finite_value(s.eigenvalue_scale);
    if (s.eigenvalue_scale <= 0 || std::cmp_not_equal(s.components.size(), n) ||
        s.component_units.size() != s.components.size() ||
        std::cmp_not_equal(s.coefficient_units.size(), square(n)))
        reject(FailureCode::InvalidInput);
    for (const double v : s.components) {
        finite_value(v);
        if (v <= 0)
            reject(FailureCode::InvalidInput);
    }
    for (core::Index j = 0; j < n; ++j)
        for (core::Index i = 0; i < n; ++i)
            for (std::size_t d = 0; d < 7; ++d) {
                const int expected =
                    static_cast<int>(s.eigenvalue_units[d]) +
                    static_cast<int>(s.component_units[static_cast<std::size_t>(i)][d]) -
                    static_cast<int>(s.component_units[static_cast<std::size_t>(j)][d]);
                if (expected !=
                    static_cast<int>(s.coefficient_units[static_cast<std::size_t>(i + j * n)][d]))
                    reject(FailureCode::InvalidInput);
            }
}
template <class T>
core::PlannedBuffer<T> buffer(core::ResourcePlan& p, core::Index n) {
    return need(p.allocate<T>(n, operation));
}
linear::BackendAdmission backend(core::Index n, const FullLimits& limits) {
    // Aggregate forecast for two Schur/Hessenberg instances, transformed/copy matrices,
    // inverse/pivot work, and temporary product packing. No stock Eigen allocation proof.
    const std::array extents{PayloadExtent{n, n, n, std::size_t{16} * 24},
                             PayloadExtent{n, 1, n, std::size_t{16} * 16}};
    const auto bytes = need(checked_payload_bytes(extents));
    return need(linear::admit_backend_memory(
        linear::detail::BackendPlanFactory::make(bytes, n == 0 ? std::optional<std::uint64_t>{0}
                                                               : std::nullopt),
        limits.backend_memory));
}
std::unique_ptr<Data> prepare(core::Index n, core::Index k, Context c, StateDiagnosticContext state,
                              const FullOutputRequest& r, core::ResourcePlan& plan, FullKind kind,
                              bool left,
                              std::optional<BasisMetricRequest> accuracy = std::nullopt) {
    auto out = std::make_unique<Data>(c, state, r, backend(k == 0 ? 0 : n, r.limits), accuracy);
    out->n = n;
    out->k = k;
    out->kind = kind;
    out->left_requested = left;
    if (r.timing == FullTimingMode::Enabled)
        out->timings.emplace();
    const auto [keep, result] = generations(r);
    out->retained_generation = keep;
    out->generation = result;
    const std::array dims{n, k};
    const auto nk = need(core::checked_extent_product(dims, operation));
    out->original = buffer<Complex>(plan, square(n));
    out->scales = buffer<double>(plan, n);
    out->units = buffer<linear::UnitDimensions>(plan, n);
    out->coefficient_units = buffer<linear::UnitDimensions>(plan, square(n));
    out->values = buffer<Complex>(plan, k);
    out->right = buffer<Complex>(plan, nk);
    if (kind == FullKind::Pairs)
        out->pairs = buffer<std::optional<FullPairProof>>(plan, k);
    else {
        out->triangular = buffer<Complex>(plan, square(k));
        out->triangular_error = buffer<double>(plan, square(k));
        out->mapping = buffer<core::Index>(plan, k);
    }
    if (left) {
        out->left = buffer<Complex>(plan, nk);
        if (kind == FullKind::Cluster) {
            out->left_triangular = buffer<Complex>(plan, square(k));
            out->left_triangular_error = buffer<double>(plan, square(k));
        }
    }
    if (std::cmp_greater(r.consumers.size(), std::numeric_limits<core::Index>::max()))
        reject(FailureCode::ResourceExceeded);
    const auto checks = static_cast<core::Index>(r.consumers.size());
    out->consumers = buffer<std::optional<core::ErrorCertificate>>(plan, checks);
    if (keep != result && kind == FullKind::Pairs)
        out->retained_pairs = buffer<std::optional<FullPairProof>>(plan, k);
    if (keep != result)
        out->retained_consumers = buffer<std::optional<core::ErrorCertificate>>(plan, checks);
    out->identity = std::make_shared<const detail::FullIdentity>();
    return out;
}
FullPacketView packet_view(const Data& d, core::NumericalGeneration gen) {
    return {d.kind,
            d.left_requested,
            d.n,
            d.k,
            d.values.values(),
            d.right.values(),
            d.left.values(),
            d.triangular.values(),
            d.left_triangular.values(),
            d.triangular_error.values(),
            d.left_triangular_error.values(),
            d.context,
            d.state,
            gen};
}
void consume(Data& d, const FullOutputRequest& r, core::NumericalGeneration gen, bool retained) {
    for (std::size_t i = 0; i < r.consumers.size(); ++i) {
        const auto check = r.consumers[i];
        core::FailureDetails details;
        details.problem = check.observation;
        details.numerical_generation = gen;
        const auto evaluated =
            core::guard_entry(operation, FailureCode::PrecisionUnresolved, details,
                              [&] { return r.evaluator->evaluate(packet_view(d, gen), check); });
        if (!evaluated) {
            auto cause = evaluated.error();
            if (!cause.details.problem)
                cause.details.problem = check.observation;
            if (!cause.details.numerical_generation)
                cause.details.numerical_generation = gen;
            throw Abort{cause};
        }
        const auto& proof = *evaluated;
        const auto& actual = proof.context();
        if (actual.observation != check.observation || actual.scope != check.scope ||
            actual.numerical_generation != gen ||
            actual.representation != core::Precision::Binary64 ||
            std::bit_cast<std::uint64_t>(proof.budget()) !=
                std::bit_cast<std::uint64_t>(check.budget) ||
            !std::holds_alternative<core::Achieved>(proof.assessment()))
            throw Abort{core::Failure{FailureCode::PrecisionUnresolved, operation, details}};
        (retained ? d.retained_consumers.values() : d.consumers.values())[i] = proof;
    }
}
struct Evaluation {
    Evaluation(Complex v, double e) noexcept : value(v), error(e) {}
    Complex value;
    double error;
};
Evaluation scaled(Complex z, double numerator, double first = 1, double second = 1) {
    const auto component = [&](double x) {
        const auto a = need(detail::magnitude(x));
        const auto b = detail::scale(a, numerator, first, second);
        if (!std::isfinite(b.upper) || (x != 0 && b.lower == 0))
            reject(FailureCode::PrecisionUnresolved);
        const double width =
            b.upper == b.lower
                ? 0.0
                : std::nextafter(b.upper - b.lower, std::numeric_limits<double>::infinity());
        return std::pair{std::copysign(b.lower, x), width};
    };
    const auto r = component(z.real()), i = component(z.imag());
    return {{r.first, i.first},
            detail::add_magnitudes({r.second, r.second}, {i.second, i.second}).upper};
}
Evaluation action_entry(const Data& d, std::span<const Complex> x, core::Index row, bool adjoint,
                        bool dimensionless) {
    detail::Sum<Complex> sum;
    double error = 0;
    for (core::Index j = 0; j < d.n; ++j) {
        const Complex a =
            adjoint ? std::conj(d.original.values()[static_cast<std::size_t>(j + row * d.n)])
                    : d.original.values()[static_cast<std::size_t>(row + j * d.n)];
        Evaluation coefficient{a, 0};
        if (dimensionless)
            coefficient = scaled(a,
                                 adjoint ? d.scales.values()[static_cast<std::size_t>(row)]
                                         : d.scales.values()[static_cast<std::size_t>(j)],
                                 adjoint ? d.scales.values()[static_cast<std::size_t>(j)]
                                         : d.scales.values()[static_cast<std::size_t>(row)],
                                 d.eigenvalue_scale);
        need(sum.product(coefficient.value, x[static_cast<std::size_t>(j)]));
        const auto magnitude = need(detail::magnitude(x[static_cast<std::size_t>(j)]));
        error = detail::add_magnitudes({error, error}, detail::multiply_magnitudes(
                                                           {coefficient.error, coefficient.error},
                                                           {magnitude.upper, magnitude.upper}))
                    .upper;
    }
    const auto value = need(sum.finish());
    error = detail::add_magnitudes({value.error, value.error}, {error, error}).upper;
    if (!finite(value.value) || !std::isfinite(error))
        reject(FailureCode::PrecisionUnresolved);
    return {value.value, error};
}
void charge(Data& d, bool adjoint) {
    auto& count = adjoint ? d.diagnostics.adjoint_applications : d.diagnostics.forward_applications;
    const auto limit = adjoint ? d.limits.actions.adjoint : d.limits.actions.forward;
    if (!count)
        count = 0;
    if (*count == std::numeric_limits<std::uint64_t>::max())
        reject(FailureCode::ResourceExceeded);
    if (*count >= limit) {
        core::FailureDetails details;
        details.resources =
            core::ResourceUsage{core::ResourceKind::OperatorApplications, *count + 1, *count};
        throw Abort{core::Failure{FailureCode::ResourceExceeded, operation, details}};
    }
    ++*count;
}
void attempt(AttemptDiagnostics& d, const FullLimits& limits) {
    if (!d.extraction_attempts)
        d.extraction_attempts = 0;
    if (*d.extraction_attempts == std::numeric_limits<std::uint64_t>::max())
        reject(FailureCode::ResourceExceeded);
    if (*d.extraction_attempts >= limits.actions.attempts) {
        core::FailureDetails details;
        details.resources = core::ResourceUsage{core::ResourceKind::Iterations,
                                                *d.extraction_attempts + 1, *d.extraction_attempts};
        throw Abort{core::Failure{FailureCode::ResourceExceeded, operation, details}};
    }
    ++*d.extraction_attempts;
}
Matrix scaled_matrix(const Data& d, bool adjoint = false) {
    Matrix a(d.n, d.n);
    for (core::Index j = 0; j < d.n; ++j)
        for (core::Index i = 0; i < d.n; ++i) {
            const auto z =
                adjoint ? std::conj(d.original.values()[static_cast<std::size_t>(j + i * d.n)])
                        : d.original.values()[static_cast<std::size_t>(i + j * d.n)];
            a(i, j) = scaled(z,
                             adjoint ? d.scales.values()[static_cast<std::size_t>(i)]
                                     : d.scales.values()[static_cast<std::size_t>(j)],
                             adjoint ? d.scales.values()[static_cast<std::size_t>(j)]
                                     : d.scales.values()[static_cast<std::size_t>(i)],
                             d.eigenvalue_scale)
                          .value;
        }
    return a;
}
void schur(const Matrix& a, Matrix& q, Matrix& t, AttemptDiagnostics& diagnostics,
           const FullLimits& limits, bool real_input) {
    if (a.rows() == 0) {
        q.resize(0, 0);
        t.resize(0, 0);
        return;
    }
    attempt(diagnostics, limits);
    if (real_input) {
        // Real 1x1 blocks have structurally real roots. Conjugate 2x2 blocks are
        // converted by a full unitary similarity, never by clamping a tiny imaginary part.
        using RealMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
        const RealMatrix real = a.real();
        Eigen::RealSchur<RealMatrix> solver(a.rows());
        solver.setMaxIterations(a.rows() * 30);
        solver.compute(real, true);
        if (solver.info() == Eigen::NoConvergence)
            reject(FailureCode::EigenCoverageUnresolved);
        if (solver.info() != Eigen::Success)
            reject(FailureCode::FactorizationFailed);
        q = solver.matrixU().cast<Complex>();
        t = solver.matrixT().cast<Complex>();
        if (!q.allFinite() || !t.allFinite())
            reject(FailureCode::FactorizationFailed);
        for (Eigen::Index i = 0; i + 1 < a.rows(); ++i) {
            if (t(i + 1, i) == Complex{})
                continue;
            const double scale =
                std::max({std::abs(t(i, i).real()), std::abs(t(i, i + 1).real()),
                          std::abs(t(i + 1, i).real()), std::abs(t(i + 1, i + 1).real())});
            if (scale == 0)
                reject(FailureCode::PrecisionUnresolved);
            const double aa = t(i, i).real() / scale, bb = t(i, i + 1).real() / scale,
                         cc = t(i + 1, i).real() / scale, dd = t(i + 1, i + 1).real() / scale;
            const double half = (aa - dd) / 2;
            const double discriminant = half * half + bb * cc;
            if (!(discriminant < 0))
                reject(FailureCode::PrecisionUnresolved);
            const Complex lambda{aa / 2 + dd / 2, std::sqrt(-discriminant)};
            Complex x = bb, y = lambda - aa;
            if (std::hypot(std::abs(lambda - dd), std::abs(cc)) >
                std::hypot(std::abs(x), std::abs(y))) {
                x = lambda - dd;
                y = cc;
            }
            Eigen::JacobiRotation<Complex> rotation;
            rotation.makeGivens(x, y);
            t.applyOnTheLeft(i, i + 1, rotation.adjoint());
            t.applyOnTheRight(i, i + 1, rotation);
            q.applyOnTheRight(i, i + 1, rotation);
            const Complex physical = lambda * scale;
            if (!finite(physical))
                reject(FailureCode::PrecisionUnresolved);
            t(i, i) = physical;
            t(i + 1, i + 1) = std::conj(physical);
            t(i + 1, i) = 0;
            ++i;
        }
    } else {
        Eigen::ComplexSchur<Matrix> solver(a.rows());
        solver.setMaxIterations(a.rows() * 30);
        solver.compute(a, true);
        if (solver.info() == Eigen::NoConvergence)
            reject(FailureCode::EigenCoverageUnresolved);
        if (solver.info() != Eigen::Success)
            reject(FailureCode::FactorizationFailed);
        q = solver.matrixU();
        t = solver.matrixT();
    }
    if (!q.allFinite() || !t.allFinite())
        reject(FailureCode::FactorizationFailed);
}
} // namespace
namespace {
WordView<Complex> words(std::span<const Complex> v) {
    return need(WordView<Complex>::make(v, {}, core::Precision::Binary64));
}
BasisMetrics basis_gate(Data& d, core::NumericalGeneration gen, bool left,
                        std::span<Complex> action, std::span<double> error) {
    if (!d.basis_request)
        reject(FailureCode::InvalidInput);
    const auto q = left ? d.left.values() : d.right.values();
    const auto t = left ? d.left_triangular.values() : d.triangular.values();
    const auto te = left ? d.left_triangular_error.values() : d.triangular_error.values();
    const auto n = static_cast<std::size_t>(d.n);
    for (core::Index j = 0; j < d.k; ++j) {
        charge(d, left);
        const auto col = q.subspan(static_cast<std::size_t>(j * d.n), n);
        for (core::Index i = 0; i < d.n; ++i) {
            const auto v = action_entry(d, col, i, left, true);
            const auto at = static_cast<std::size_t>(i + j * d.n);
            action[at] = v.value;
            error[at] = v.error;
        }
    }
    return need(assess_schur_metrics(d.n, d.k, words(q), ActionValues<Complex>{words(t), te},
                                     ActionValues<Complex>{words(action), error}, *d.basis_request,
                                     MetricContext{d.context, d.state, gen}));
}
void qualify_basis(Data& d, const FullRequest& request, std::span<Complex> action,
                   std::span<double> error) {
    d.basis = basis_gate(d, d.retained_generation, false, action, error);
    if (d.left_requested)
        d.left_basis = basis_gate(d, d.retained_generation, true, action, error);
    const bool projected = d.generation != d.retained_generation;
    if (projected) {
        d.retained_basis = d.basis;
        d.retained_left_basis = d.left_basis;
    }
    consume(d, request, d.retained_generation, projected);
    if (projected) {
        d.basis = basis_gate(d, d.generation, false, action, error);
        if (d.left_requested)
            d.left_basis = basis_gate(d, d.generation, true, action, error);
        consume(d, request, d.generation, false);
    }
}
void copy_snapshot(Data& out, const Data& source) {
    std::copy(source.original.values().begin(), source.original.values().end(),
              out.original.values().begin());
    std::copy(source.scales.values().begin(), source.scales.values().end(),
              out.scales.values().begin());
    std::copy(source.units.values().begin(), source.units.values().end(),
              out.units.values().begin());
    std::copy(source.coefficient_units.values().begin(), source.coefficient_units.values().end(),
              out.coefficient_units.values().begin());
    out.eigenvalue_scale = source.eigenvalue_scale;
    out.eigenvalue_units = source.eigenvalue_units;
    out.real_input = source.real_input;
}
void import_basis(Data& d, const Matrix& q, const Matrix& t, bool left) {
    auto target = left ? d.left.values() : d.right.values();
    auto triangular = left ? d.left_triangular.values() : d.triangular.values();
    for (core::Index j = 0; j < d.k; ++j) {
        for (core::Index i = 0; i < d.n; ++i)
            target[static_cast<std::size_t>(i + j * d.n)] = q(i, j);
        for (core::Index i = 0; i < d.k; ++i)
            triangular[static_cast<std::size_t>(i + j * d.k)] = i <= j ? t(i, j) : Complex{};
    }
}
void physical_diagonal(Data& d, bool left) {
    auto t = left ? d.left_triangular.values() : d.triangular.values();
    auto errors = left ? d.left_triangular_error.values() : d.triangular_error.values();
    for (core::Index j = 0; j < d.k; ++j) {
        const auto sorted = std::find(d.mapping.values().begin(), d.mapping.values().end(), j) -
                            d.mapping.values().begin();
        if (std::cmp_greater_equal(sorted, d.values.size()))
            reject(FailureCode::InvalidInput);
        const auto lambda = left ? std::conj(d.values.values()[static_cast<std::size_t>(sorted)])
                                 : d.values.values()[static_cast<std::size_t>(sorted)];
        const auto actual = scaled(lambda, 1, d.eigenvalue_scale);
        const auto at = static_cast<std::size_t>(j + j * d.k);
        t[at] = actual.value;
        errors[at] = actual.error;
    }
}
void initialize_diagnostics(AttemptDiagnostics& d, core::Index count) {
    d.requested_count = static_cast<std::uint64_t>(count);
    d.extraction_attempts = 0;
    d.forward_applications = 0;
    d.adjoint_applications = 0;
    d.inverse_applications = 0;
}
struct ParentGuard {
    const Data& data;
    explicit ParentGuard(const Data& d) : data(d) {
        bool expected = false;
        if (!data.busy.compare_exchange_strong(expected, true))
            reject(FailureCode::InvalidInput);
    }
    ~ParentGuard() {
        data.busy.store(false);
    }
};
} // namespace

namespace detail {
template <StorageScalar Scalar>
eigen::Result<FullSpectrum> FullFactory::spectrum(const PrecisionBuffer<Scalar>& original,
                                                  Context context, StateDiagnosticContext state,
                                                  FullScales scales, const FullRequest& request,
                                                  core::ResourcePlan& plan) noexcept {
    AttemptDiagnostics diagnostics;
    return boundary<FullSpectrum>(context, state, diagnostics, [&]() -> FullSpectrum {
        const FullPreparationClock preparation{request.timing};
        const auto shape = original.shape();
        if (shape.rows < 0 || shape.rows != shape.columns)
            reject(FailureCode::InvalidInput);
        const auto n = shape.rows;
        initialize_diagnostics(diagnostics, n);
        if (original.precision() != core::Precision::Binary64)
            reject(FailureCode::PrecisionUnsupported);
        if (original.generation() != context.cache.matrix)
            reject(FailureCode::StaleGeneration);
        const auto* storage = std::get_if<Binary64<Scalar>>(&original.storage());
        if (!storage)
            reject(FailureCode::PrecisionUnsupported);
        need(storage->values.validate(operation));
        preflight(n, context, request, n);
        scale_inputs(n, scales);
        const auto [keep, returned] = generations(request);
        if (keep == original.generation() || returned == original.generation())
            reject(FailureCode::InvalidInput);
        auto data =
            prepare(n, n, context, state, request, plan, FullKind::Spectrum, false, request.basis);
        data->diagnostics = diagnostics;
        const ReportSync report{diagnostics, *data};
        data->real_input = std::same_as<Scalar, double>;
        data->eigenvalue_scale = scales.eigenvalue_scale;
        data->eigenvalue_units = scales.eigenvalue_units;
        std::copy(scales.components.begin(), scales.components.end(),
                  data->scales.values().begin());
        std::copy(scales.component_units.begin(), scales.component_units.end(),
                  data->units.values().begin());
        std::copy(scales.coefficient_units.begin(), scales.coefficient_units.end(),
                  data->coefficient_units.values().begin());
        for (core::Index j = 0; j < n; ++j)
            for (core::Index i = 0; i < n; ++i)
                data->original.values()[static_cast<std::size_t>(i + j * n)] =
                    storage->values
                        .values()[static_cast<std::size_t>(i + j * shape.leading_dimension)];
        auto actions = buffer<Complex>(plan, square(n));
        auto errors = buffer<double>(plan, square(n));
        const Matrix a = scaled_matrix(*data);
        Matrix q(n, n), t(n, n);
        Matrix physical(n, 1);
        preparation.finish(*data);
        if (n > 0) {
            const FullStageClock timing{*data, &FullTimings::backend};
            schur(a, q, t, data->diagnostics, request.limits, data->real_input);
        }
        {
            const FullStageClock qualification{*data, &FullTimings::qualification};
            data->diagnostics.converged_count = static_cast<std::uint64_t>(n);
            import_basis(*data, q, t, false);
            // The backend factors are candidates only. The final gate below evaluates
            // the actual physical eigenvalues and the corresponding returned T diagonal.
            for (core::Index i = 0; i < n; ++i) {
                physical(i) = scaled(t(i, i), scales.eigenvalue_scale).value;
                data->mapping.values()[static_cast<std::size_t>(i)] = i;
            }
            std::iota(data->mapping.values().begin(), data->mapping.values().end(), core::Index{0});
            std::sort(data->mapping.values().begin(), data->mapping.values().end(),
                      [&](core::Index i, core::Index j) {
                          const auto x = physical(i), y = physical(j);
                          return x.real() != y.real() ? x.real() < y.real() : x.imag() < y.imag();
                      });
            for (core::Index i = 0; i < n; ++i)
                data->values.values()[static_cast<std::size_t>(i)] =
                    physical(data->mapping.values()[static_cast<std::size_t>(i)]);
            physical_diagonal(*data, false);
            qualify_basis(*data, request, actions.values(), errors.values());
        }
        check_timing(*data);
        return FullSpectrum{std::move(data)};
    });
}
eigen::Result<FullSelection> FullFactory::selection(const FullSpectrum& packet,
                                                    std::span<const core::Index> indices,
                                                    core::ResourcePlan& plan) noexcept {
    if (!packet.data_)
        return std::unexpected(
            EigenFailure{core::Failure{FailureCode::InvalidInput, operation}, {}});
    const Context context = packet.data_->context;
    const StateDiagnosticContext state = packet.data_->state;
    AttemptDiagnostics diagnostics;
    return boundary<FullSelection>(context, state, diagnostics, [&]() -> FullSelection {
        if (!packet.data_)
            reject(FailureCode::InvalidInput);
        const auto& d = *packet.data_;
        const ParentGuard guard{d};
        if (std::cmp_greater(indices.size(), d.n))
            reject(FailureCode::InvalidInput);
        for (std::size_t i = 0; i < indices.size(); ++i) {
            if (indices[i] < 0 || indices[i] >= d.n)
                reject(FailureCode::InvalidInput);
            for (std::size_t j = 0; j < i; ++j)
                if (indices[i] == indices[j])
                    reject(FailureCode::InvalidInput);
        }
        auto copy = buffer<core::Index>(plan, static_cast<core::Index>(indices.size()));
        std::copy(indices.begin(), indices.end(), copy.values().begin());
        std::sort(copy.values().begin(), copy.values().end());
        return FullSelection{d.identity, std::move(copy)};
    });
}
} // namespace detail
FullSelection::FullSelection(std::shared_ptr<const detail::FullIdentity> identity,
                             core::PlannedBuffer<core::Index> indices) noexcept
    : identity_(std::move(identity)), indices_(std::move(indices)) {}
Result<FullSelection> FullSelection::make(const FullSpectrum& p, std::span<const core::Index> i,
                                          core::ResourcePlan& plan) noexcept {
    return detail::FullFactory::selection(p, i, plan);
}
template <class Q>
FullPacket<Q>::FullPacket(std::unique_ptr<Data> data) noexcept : data_(std::move(data)) {}
template <class Q>
FullPacket<Q>::FullPacket(FullPacket&&) noexcept = default;
template <class Q>
FullPacket<Q>& FullPacket<Q>::operator=(FullPacket&&) noexcept = default;
template <class Q>
FullPacket<Q>::~FullPacket() = default;
template <class Q>
bool FullPacket<Q>::active() const noexcept {
    return data_ != nullptr;
}
template <class Q>
std::optional<FullPacketView> FullPacket<Q>::view() const& noexcept {
    if (!data_)
        return std::nullopt;
    return packet_view(*data_, data_->generation);
}
template <class Q>
AttemptDiagnostics FullPacket<Q>::diagnostics() const noexcept {
    return data_ ? data_->diagnostics : AttemptDiagnostics{};
}
template <class Q>
FullPhysicalCoverage FullPacket<Q>::physical_coverage() const noexcept {
    return FullPhysicalCoverage::Unresolved;
}
template <class Q>
const linear::ReturnedPrecision* FullPacket<Q>::returned_precision() const& noexcept {
    return data_ ? &data_->returned : nullptr;
}
template <class Q>
std::span<const core::Index> FullPacket<Q>::diagonal_mapping() const& noexcept {
    return data_ ? data_->mapping.values() : std::span<const core::Index>{};
}
template <class Q>
const BasisMetrics* FullPacket<Q>::basis_proof() const& noexcept {
    return data_ && data_->basis ? &*data_->basis : nullptr;
}
template <class Q>
const BasisMetrics* FullPacket<Q>::left_basis_proof() const& noexcept {
    return data_ && data_->left_basis ? &*data_->left_basis : nullptr;
}
template <class Q>
Result<FullPairProof> FullPacket<Q>::pair_proof(core::Index index) const noexcept {
    if (!data_)
        return std::unexpected(
            EigenFailure{core::Failure{FailureCode::InvalidInput, operation}, {}});
    if (index < 0 || std::cmp_greater_equal(index, data_->pairs.size()))
        return std::unexpected(complete_failure(core::Failure{FailureCode::InvalidInput, operation},
                                                data_->context, data_->state));
    const auto& proof = data_->pairs.values()[static_cast<std::size_t>(index)];
    if (!proof)
        return std::unexpected(
            complete_failure(core::Failure{FailureCode::PrecisionUnresolved, operation},
                             data_->context, data_->state, data_->diagnostics));
    return *proof;
}
template <class Q>
std::size_t FullPacket<Q>::consumer_count() const noexcept {
    return data_ ? data_->consumers.size() : 0;
}
template <class Q>
Result<core::ErrorCertificate> FullPacket<Q>::consumer_proof(core::Index index) const noexcept {
    if (!data_)
        return std::unexpected(
            EigenFailure{core::Failure{FailureCode::InvalidInput, operation}, {}});
    if (index < 0 || std::cmp_greater_equal(index, data_->consumers.size()))
        return std::unexpected(complete_failure(core::Failure{FailureCode::InvalidInput, operation},
                                                data_->context, data_->state));
    const auto& proof = data_->consumers.values()[static_cast<std::size_t>(index)];
    if (!proof)
        return std::unexpected(
            complete_failure(core::Failure{FailureCode::PrecisionUnresolved, operation},
                             data_->context, data_->state, data_->diagnostics));
    return *proof;
}
template <StorageScalar S>
Result<FullSpectrum> full_spectrum(const PrecisionBuffer<S>& a, Context c, StateDiagnosticContext s,
                                   FullScales scales, const FullRequest& r,
                                   core::ResourcePlan& p) noexcept {
    return detail::FullFactory::spectrum(a, c, s, scales, r, p);
}
template <class Q>
std::optional<FullTimings> FullPacket<Q>::timings() const noexcept {
    return data_ ? data_->timings : std::nullopt;
}
template <class Q>
std::uint64_t FullPacket<Q>::clock_reads() const noexcept {
    return data_ ? data_->clock_reads : 0;
}
template <class Q>
std::optional<double> FullPacket<Q>::independence_upper() const noexcept {
    return data_ ? data_->independence_upper : std::nullopt;
}
template <class Q>
std::optional<FullScales> FullPacket<Q>::scales() const& noexcept {
    if (!data_)
        return std::nullopt;
    return FullScales{data_->scales.values(), data_->eigenvalue_scale, data_->units.values(),
                      data_->eigenvalue_units, data_->coefficient_units.values()};
}
template <class Q>
std::span<const Complex> FullPacket<Q>::original_coefficients() const& noexcept {
    return data_ ? data_->original.values() : std::span<const Complex>{};
}
template <class Q>
std::optional<FullLimits> FullPacket<Q>::limits() const noexcept {
    return data_ ? std::optional<FullLimits>{data_->limits} : std::nullopt;
}
template <class Q>
std::optional<linear::BackendAdmission> FullPacket<Q>::backend_admission() const noexcept {
    return data_ ? std::optional<linear::BackendAdmission>{data_->backend} : std::nullopt;
}
template <class Q>
std::optional<FullAssociation> FullPacket<Q>::association_policy() const noexcept {
    return data_ ? data_->association : std::nullopt;
}
template <class Q>
std::optional<FullBackend> FullPacket<Q>::source_backend() const noexcept {
    if (!data_)
        return std::nullopt;
    return data_->real_input ? FullBackend::RealSchurComplexified : FullBackend::ComplexSchur;
}

template <class Q>
std::optional<core::Index> FullPacket<Q>::schur_iteration_limit() const noexcept {
    return data_ && data_->kind == FullKind::Spectrum ? std::optional<core::Index>{30 * data_->n}
                                                      : std::nullopt;
}
template class FullPacket<FullSpectrumQualified>;
template class FullPacket<FullPairsQualified>;
template class FullPacket<FullClusterQualified>;
template Result<FullSpectrum> full_spectrum(const PrecisionBuffer<double>&, Context,
                                            StateDiagnosticContext, FullScales, const FullRequest&,
                                            core::ResourcePlan&) noexcept;
template Result<FullSpectrum> full_spectrum(const PrecisionBuffer<Complex>&, Context,
                                            StateDiagnosticContext, FullScales, const FullRequest&,
                                            core::ResourcePlan&) noexcept;
} // namespace uvlas::numerics::eigen
namespace uvlas::numerics::eigen {
namespace {
void reorder(Matrix& q, Matrix& t, std::span<const core::Index> desired,
             std::span<core::Index> positions) {
    std::iota(positions.begin(), positions.end(), core::Index{0});
    for (std::size_t target = 0; target < desired.size(); ++target) {
        auto found = std::find(positions.begin() + static_cast<std::ptrdiff_t>(target),
                               positions.end(), desired[target]);
        if (found == positions.end())
            reject(FailureCode::InvalidInput);
        auto index = static_cast<core::Index>(found - positions.begin());
        while (index > static_cast<core::Index>(target)) {
            const auto p = index - 1;
            const Complex a = t(p, p), b = t(p, p + 1), d = t(p + 1, p + 1);
            const double scale =
                std::max({std::abs(a.real()), std::abs(a.imag()), std::abs(b.real()),
                          std::abs(b.imag()), std::abs(d.real()), std::abs(d.imag())});
            Eigen::JacobiRotation<Complex> rotation;
            if (scale > 0 && (b != Complex{} || a != d))
                rotation.makeGivens(b / scale, d / scale - a / scale);
            else
                rotation = Eigen::JacobiRotation<Complex>{0, 1};
            // Apply the full similarity, including couplings outside the adjacent block.
            t.applyOnTheLeft(p, p + 1, rotation.adjoint());
            t.applyOnTheRight(p, p + 1, rotation);
            q.applyOnTheRight(p, p + 1, rotation);
            if (!q.allFinite() || !t.allFinite())
                reject(FailureCode::PrecisionUnresolved);
            // This entry is exactly zero under the unitary exchange. Its discarded
            // floating residue is never used as a qualification: original A is rechecked.
            t(p + 1, p) = 0;
            std::swap(positions[static_cast<std::size_t>(p)],
                      positions[static_cast<std::size_t>(p + 1)]);
            --index;
        }
    }
}
bool associated(Complex x, Complex y, FullAssociation policy) {
    detail::Sum<Complex> difference;
    need(difference.add(x));
    need(difference.add(y, -1));
    const auto delta = need(difference.finish());
    const auto distance = need(detail::magnitude(delta.value));
    const double upper = detail::add_magnitudes(distance, {delta.error, delta.error}).upper;
    const double lower = detail::lower_subtract(distance.lower, delta.error);
    const auto ax = need(detail::magnitude(x)), ay = need(detail::magnitude(y));
    const auto size = detail::Bounds{std::max(ax.lower, ay.lower), std::max(ax.upper, ay.upper)};
    const auto allowance = detail::add_magnitudes(
        {policy.absolute, policy.absolute},
        detail::multiply_magnitudes({policy.relative, policy.relative}, size));
    if (upper <= allowance.lower)
        return true;
    if (lower > allowance.upper)
        return false;
    reject(FailureCode::PrecisionUnresolved);
}
void adjoint_schur(const Data& parent, Matrix& q, Matrix& t) {
    const auto n = parent.n;
    q.resize(n, n);
    t.resize(n, n);
    // Ahat*=Q T* Q*. Reversing both Schur coordinates turns T* into an
    // upper triangular factor without perturbing a defective block's eigenvalues.
    // A subsequent left-specific exchange and original A* qualification are mandatory.
    for (core::Index j = 0; j < n; ++j) {
        for (core::Index i = 0; i < n; ++i) {
            q(i, j) = parent.right.values()[static_cast<std::size_t>(i + (n - 1 - j) * n)];
            t(i, j) =
                std::conj(parent.triangular
                              .values()[static_cast<std::size_t>((n - 1 - j) + (n - 1 - i) * n)]);
        }
    }
}
void diagonal_matches(const Data& source, const Matrix& t, std::span<const core::Index> selected,
                      FullAssociation policy, std::span<core::Index> matches,
                      std::span<core::Index> owner, std::span<core::Index> visited) {
    std::fill(owner.begin(), owner.end(), core::Index{-1});
    const auto augment = [&](auto&& self, core::Index item, core::Index stamp) -> bool {
        const auto lambda = scaled(source.values.values()[static_cast<std::size_t>(
                                       selected[static_cast<std::size_t>(item)])],
                                   1, source.eigenvalue_scale)
                                .value;
        for (core::Index j = 0; j < source.n; ++j) {
            if (visited[static_cast<std::size_t>(j)] == stamp ||
                !associated(std::conj(lambda), t(j, j), policy))
                continue;
            visited[static_cast<std::size_t>(j)] = stamp;
            const auto previous = owner[static_cast<std::size_t>(j)];
            if (previous < 0 || self(self, previous, stamp)) {
                owner[static_cast<std::size_t>(j)] = item;
                matches[static_cast<std::size_t>(item)] = j;
                return true;
            }
        }
        return false;
    };
    // Match whole clusters with multiplicity, without assigning identities inside repeated roots.
    std::fill(visited.begin(), visited.end(), core::Index{-1});
    for (core::Index i = 0; std::cmp_less(i, selected.size()); ++i)
        if (!augment(augment, i, i))
            reject(FailureCode::PrecisionUnresolved);
}
void child_generations(const Data& parent, const FullOutputRequest& r) {
    const auto [retained, returned] = generations(r);
    if (retained == parent.generation || retained == parent.retained_generation ||
        returned == parent.generation || returned == parent.retained_generation)
        reject(FailureCode::InvalidInput);
}
} // namespace
namespace detail {
eigen::Result<FullCluster> FullFactory::cluster(const FullSpectrum& packet,
                                                const FullSelection& selection,
                                                const FullClusterRequest& request,
                                                core::ResourcePlan& plan) noexcept {
    if (!packet.data_)
        return std::unexpected(
            EigenFailure{core::Failure{FailureCode::InvalidInput, operation}, {}});
    const auto& parent = *packet.data_;
    AttemptDiagnostics diagnostics;
    return boundary<FullCluster>(parent.context, parent.state, diagnostics, [&]() -> FullCluster {
        const FullPreparationClock preparation{request.common.timing};
        const ParentGuard guard{parent};
        if (selection.identity_ != parent.identity)
            reject(FailureCode::StaleGeneration);
        const auto n = parent.n, k = static_cast<core::Index>(selection.indices_.size());
        initialize_diagnostics(diagnostics, k);
        preflight(n, parent.context, request.common, k);
        association(request.association);
        child_generations(parent, request.common);
        auto data = prepare(n, k, parent.context, parent.state, request.common, plan,
                            FullKind::Cluster, request.left, request.common.basis);
        data->diagnostics = diagnostics;
        const ReportSync report{diagnostics, *data};
        copy_snapshot(*data, parent);
        data->association = request.association;
        const std::array dims{n, k};
        const auto count = need(core::checked_extent_product(dims, operation));
        auto actions = buffer<Complex>(plan, count);
        auto errors = buffer<double>(plan, count);
        auto desired = buffer<core::Index>(plan, k);
        auto positions = buffer<core::Index>(plan, n);
        auto matched = buffer<core::Index>(plan, k);
        auto owner = buffer<core::Index>(plan, n);
        auto visited = buffer<core::Index>(plan, n);
        for (core::Index j = 0; j < k; ++j) {
            const auto at = static_cast<std::size_t>(j),
                       selected = static_cast<std::size_t>(selection.indices_.values()[at]);
            desired.values()[at] = parent.mapping.values()[selected];
            data->values.values()[at] = parent.values.values()[selected];
            data->mapping.values()[at] = j;
        }
        preparation.finish(*data);
        if (k > 0) {
            const FullStageClock timing{*data, &FullTimings::backend};
            attempt(data->diagnostics, request.common.limits);
            Matrix q(n, n), t(n, n);
            std::copy(parent.right.values().begin(), parent.right.values().end(), q.data());
            std::copy(parent.triangular.values().begin(), parent.triangular.values().end(),
                      t.data());
            // Parent diagonal includes its physical-value conversion center. Swaps produce
            // candidates only; unchanged selected physical lambda values are authoritative.
            reorder(q, t, desired.values(), positions.values());
            import_basis(*data, q, t, false);
            physical_diagonal(*data, false);
            if (request.left) {
                Matrix leftq(n, n), leftt(n, n);
                adjoint_schur(parent, leftq, leftt);
                diagonal_matches(parent, leftt, selection.indices_.values(), request.association,
                                 matched.values(), owner.values(), visited.values());
                reorder(leftq, leftt, matched.values(), positions.values());
                import_basis(*data, leftq, leftt, true);
                physical_diagonal(*data, true);
            }
        }
        {
            const FullStageClock qualification{*data, &FullTimings::qualification};
            data->diagnostics.converged_count = static_cast<std::uint64_t>(k);
            qualify_basis(*data, request.common, actions.values(), errors.values());
        }
        check_timing(*data);
        return FullCluster{std::move(data)};
    });
}
} // namespace detail
Result<FullCluster> full_cluster(const FullSpectrum& packet, const FullSelection& selected,
                                 const FullClusterRequest& request,
                                 core::ResourcePlan& plan) noexcept {
    return detail::FullFactory::cluster(packet, selected, request, plan);
}
} // namespace uvlas::numerics::eigen
namespace uvlas::numerics::eigen {
namespace {
void eigenvector(const Matrix& q, const Matrix& t, core::Index index, std::span<Complex> work,
                 std::span<Complex> result) {
    std::fill(work.begin(), work.end(), Complex{});
    work[static_cast<std::size_t>(index)] = 1;
    for (core::Index row = index; row-- > 0;) {
        detail::Sum<Complex> sum;
        for (core::Index j = row + 1; j <= index; ++j)
            need(sum.product(t(row, j), work[static_cast<std::size_t>(j)]));
        const auto value = need(sum.finish());
        const Complex denominator = t(row, row) - t(index, index);
        // No rank tolerance or artificial diagonal shift is introduced. An unresolved
        // triangular direction is a failed individual-pair request, not a failed spectrum.
        if (denominator == Complex{})
            reject(FailureCode::PrecisionUnresolved);
        const auto coefficient = -value.value / denominator;
        if (!finite(coefficient))
            reject(FailureCode::PrecisionUnresolved);
        work[static_cast<std::size_t>(row)] = coefficient;
    }
    double largest = 0;
    for (const auto value : work)
        largest = std::max(largest, std::max(std::abs(value.real()), std::abs(value.imag())));
    if (largest == 0 || !std::isfinite(largest))
        reject(FailureCode::PrecisionUnresolved);
    for (auto& value : work)
        value /= largest;
    for (core::Index i = 0; i < q.rows(); ++i) {
        detail::Sum<Complex> sum;
        for (core::Index j = 0; j < q.cols(); ++j)
            need(sum.product(q(i, j), work[static_cast<std::size_t>(j)]));
        result[static_cast<std::size_t>(i)] = need(sum.finish()).value;
    }
    double norm = 0;
    std::size_t pivot = 0;
    double peak = 0;
    for (std::size_t i = 0; i < result.size(); ++i) {
        norm = std::hypot(norm, std::abs(result[i]));
        const double a = std::abs(result[i]);
        if (a > peak) {
            peak = a;
            pivot = i;
        }
    }
    if (norm == 0 || !std::isfinite(norm))
        reject(FailureCode::PrecisionUnresolved);
    const Complex phase = std::conj(result[pivot]) / std::abs(result[pivot]);
    for (auto& value : result)
        value = (value / norm) * phase;
}
void physical_phase(std::span<Complex> physical, std::span<const double> scales, PairSide side,
                    std::span<Complex> backup, std::span<Complex> candidate) {
    std::copy(physical.begin(), physical.end(), backup.begin());
    // Rotation itself rounds. Test each nonzero phase anchor against the exact
    // maximum of its actual returned physical words, rather than assuming a pivot
    // selected before D conversion still wins afterward. No fuzzy tie policy is used.
    for (std::size_t pivot = 0; pivot < physical.size(); ++pivot) {
        const auto z = backup[pivot];
        const double largest = std::max(std::abs(z.real()), std::abs(z.imag()));
        if (largest == 0)
            continue;
        const auto normalized = z / largest;
        const auto phase = std::conj(normalized) / std::abs(normalized);
        const double magnitude = std::abs(z);
        if (!std::isfinite(magnitude))
            continue;
        bool valid = true;
        for (std::size_t i = 0; i < physical.size(); ++i) {
            candidate[i] = backup[i] * phase;
            if (!finite(candidate[i]))
                valid = false;
        }
        if (!valid)
            continue;
        candidate[pivot] = {magnitude, 0};
        std::size_t maximum = 0;
        for (std::size_t i = 1; i < physical.size(); ++i)
            if (detail::full_phase_coordinate_larger(candidate[i], scales[i], candidate[maximum],
                                                     scales[maximum], side))
                maximum = i;
        if (maximum == pivot) {
            std::copy(candidate.begin(), candidate.end(), physical.begin());
            return;
        }
    }
    reject(FailureCode::PrecisionUnresolved);
}
void isolated(const Data& parent, core::Index selected, FullAssociation policy, bool full_basis) {
    const auto lambda = scaled(parent.values.values()[static_cast<std::size_t>(selected)], 1,
                               parent.eigenvalue_scale)
                            .value;
    for (core::Index other = 0; other < parent.n; ++other) {
        if (other == selected)
            continue;
        const auto other_value = scaled(parent.values.values()[static_cast<std::size_t>(other)], 1,
                                        parent.eigenvalue_scale)
                                     .value;
        if (associated(lambda, other_value, policy))
            reject(full_basis ? FailureCode::RankUnresolved : FailureCode::PrecisionUnresolved);
    }
}
PairMetrics pair_gate(Data& d, core::Index column, bool left, core::NumericalGeneration generation,
                      const PairMetricRequest& request, std::span<Complex> action,
                      std::span<double> error) {
    const auto start = static_cast<std::size_t>(column * d.n), n = static_cast<std::size_t>(d.n);
    const auto vector = (left ? d.left.values() : d.right.values()).subspan(start, n);
    charge(d, left);
    for (core::Index i = 0; i < d.n; ++i) {
        const auto value = action_entry(d, vector, i, left, false);
        action[static_cast<std::size_t>(i)] = value.value;
        error[static_cast<std::size_t>(i)] = value.error;
    }
    return need(assess_pair_metrics(
        words(vector), words(d.values.values().subspan(static_cast<std::size_t>(column), 1)),
        ActionValues<Complex>{words(action), error},
        MetricScales{d.scales.values(), d.eigenvalue_scale},
        left ? PairSide::Left : PairSide::Right, request,
        MetricContext{d.context, d.state, generation}));
}
FullSensitivity sensitivity(const Data& d, core::Index column) {
    detail::Sum<Complex> dot;
    double error = 0, right_norm = 0, left_norm = 0;
    for (core::Index i = 0; i < d.n; ++i) {
        const auto at = static_cast<std::size_t>(i + column * d.n);
        const auto z =
            scaled(d.right.values()[at], 1, d.scales.values()[static_cast<std::size_t>(i)]);
        const auto y = scaled(d.left.values()[at], d.scales.values()[static_cast<std::size_t>(i)]);
        need(dot.product(std::conj(y.value), z.value));
        const auto az = need(detail::magnitude(z.value)), ay = need(detail::magnitude(y.value));
        const auto first = detail::multiply_magnitudes({y.error, y.error}, az);
        const auto second = detail::multiply_magnitudes({z.error, z.error}, ay);
        const auto cross = detail::multiply_magnitudes({y.error, y.error}, {z.error, z.error});
        error = detail::add_magnitudes(
                    {error, error},
                    detail::add_magnitudes(first, detail::add_magnitudes(second, cross)))
                    .upper;
        right_norm = std::hypot(right_norm, std::abs(z.value));
        left_norm = std::hypot(left_norm, std::abs(y.value));
    }
    const auto value = need(dot.finish());
    const auto mag = need(detail::magnitude(value.value));
    error = detail::add_magnitudes({error, error}, {value.error, value.error}).upper;
    if (detail::lower_subtract(mag.lower, error) <= 0)
        return {std::nullopt, true};
    const double estimate = (right_norm * left_norm) / std::abs(value.value);
    if (!std::isfinite(estimate))
        return {std::nullopt, true};
    return {estimate, false};
}
void independence(Data& d, Matrix& z, Matrix& inverse) {
    if (d.n == 0) {
        d.independence_upper = 0;
        return;
    }
    for (core::Index j = 0; j < d.n; ++j)
        for (core::Index i = 0; i < d.n; ++i)
            z(i, j) = scaled(d.right.values()[static_cast<std::size_t>(i + j * d.n)], 1,
                             d.scales.values()[static_cast<std::size_t>(i)])
                          .value;
    const Eigen::PartialPivLU<Matrix> lu(z);
    const auto diag = lu.matrixLU().diagonal();
    for (core::Index i = 0; i < d.n; ++i)
        if (diag[i] == Complex{} || !finite(diag[i]))
            reject(FailureCode::RankUnresolved);
    inverse = lu.inverse();
    if (!inverse.allFinite())
        reject(FailureCode::RankUnresolved);
    double maximum = 0;
    for (core::Index i = 0; i < d.n; ++i) {
        double row_upper = 0;
        for (core::Index j = 0; j < d.n; ++j) {
            detail::Sum<Complex> sum;
            need(sum.add(i == j ? Complex{1} : Complex{}));
            double error = 0;
            for (core::Index k = 0; k < d.n; ++k) {
                const auto actual = scaled(d.right.values()[static_cast<std::size_t>(k + j * d.n)],
                                           1, d.scales.values()[static_cast<std::size_t>(k)]);
                need(sum.product(inverse(i, k), actual.value, -1));
                const auto a = need(detail::magnitude(inverse(i, k)));
                error = detail::add_magnitudes({error, error}, detail::multiply_magnitudes(
                                                                   a, {actual.error, actual.error}))
                            .upper;
            }
            const auto value = need(sum.finish());
            const auto absolute = need(detail::magnitude(value.value));
            const auto upper =
                detail::add_magnitudes(
                    absolute, detail::add_magnitudes({value.error, value.error}, {error, error}))
                    .upper;
            row_upper = detail::add_magnitudes({row_upper, row_upper}, {upper, upper}).upper;
        }
        maximum = std::max(maximum, row_upper);
    }
    if (!std::isfinite(maximum) || maximum >= 1)
        reject(FailureCode::RankUnresolved);
    d.independence_upper = maximum;
}
void qualify_pairs(Data& d, const FullPairRequest& r, std::span<Complex> action,
                   std::span<double> error) {
    const auto phase = [&](core::NumericalGeneration generation, bool retained) {
        for (core::Index j = 0; j < d.k; ++j) {
            FullPairProof proof{pair_gate(d, j, false, generation, r.pair, action, error),
                                std::nullopt, std::nullopt};
            if (r.left) {
                proof.left = pair_gate(d, j, true, generation, r.pair, action, error);
                proof.sensitivity = sensitivity(d, j);
            }
            (retained ? d.retained_pairs.values() : d.pairs.values())[static_cast<std::size_t>(j)] =
                proof;
        }
        consume(d, r.common, generation, retained);
    };
    const bool projection = d.retained_generation != d.generation;
    phase(d.retained_generation, projection);
    if (projection)
        phase(d.generation, false);
}
} // namespace
namespace detail {
eigen::Result<FullPairs> FullFactory::pairs(const FullSpectrum& packet,
                                            const FullSelection& selection,
                                            const FullPairRequest& request,
                                            core::ResourcePlan& plan) noexcept {
    if (!packet.data_)
        return std::unexpected(
            EigenFailure{core::Failure{FailureCode::InvalidInput, operation}, {}});
    const auto& parent = *packet.data_;
    AttemptDiagnostics diagnostics;
    return boundary<FullPairs>(parent.context, parent.state, diagnostics, [&]() -> FullPairs {
        const FullPreparationClock preparation{request.common.timing};
        const ParentGuard guard{parent};
        if (selection.identity_ != parent.identity)
            reject(FailureCode::StaleGeneration);
        const auto n = parent.n, k = static_cast<core::Index>(selection.indices_.size());
        initialize_diagnostics(diagnostics, k);
        output_preflight(n, parent.context, request.common, k);
        association(request.association);
        child_generations(parent, request.common);
        requirement(request.pair.residual);
        requirement(request.pair.normalization);
        backward_requirement(request.pair.backward);
        if (request.pair.residual.budget > 1e-10)
            reject(FailureCode::InvalidInput);
        if (request.require_full_basis && k != n)
            reject(FailureCode::InvalidInput);
        for (auto index : selection.indices_.values())
            isolated(parent, index, request.association, request.require_full_basis);
        auto data = prepare(n, k, parent.context, parent.state, request.common, plan,
                            FullKind::Pairs, request.left);
        data->diagnostics = diagnostics;
        const ReportSync report{diagnostics, *data};
        copy_snapshot(*data, parent);
        data->association = request.association;
        auto work = buffer<Complex>(plan, n);
        auto candidate = buffer<Complex>(plan, n);
        auto actions = buffer<Complex>(plan, n);
        auto errors = buffer<double>(plan, n);
        const auto work_order = k > 0 ? n : 0;
        Matrix q(work_order, work_order), t(work_order, work_order),
            leftq(request.left ? work_order : 0, request.left ? work_order : 0),
            leftt(request.left ? work_order : 0, request.left ? work_order : 0);
        Matrix inverse(request.require_full_basis ? n : 0, request.require_full_basis ? n : 0),
            z(inverse.rows(), inverse.cols());
        preparation.finish(*data);
        if (k > 0) {
            const FullStageClock timing{*data, &FullTimings::backend};
            attempt(data->diagnostics, request.common.limits);
            std::copy(parent.right.values().begin(), parent.right.values().end(), q.data());
            std::copy(parent.triangular.values().begin(), parent.triangular.values().end(),
                      t.data());
            if (request.left)
                adjoint_schur(parent, leftq, leftt);
            for (core::Index j = 0; j < k; ++j) {
                const auto at = static_cast<std::size_t>(j),
                           selected = static_cast<std::size_t>(selection.indices_.values()[at]);
                const auto diagonal = parent.mapping.values()[selected];
                data->values.values()[at] = parent.values.values()[selected];
                eigenvector(q, t, diagonal, work.values(), candidate.values());
                for (core::Index i = 0; i < n; ++i)
                    data->right.values()[static_cast<std::size_t>(i + j * n)] =
                        scaled(candidate.values()[static_cast<std::size_t>(i)],
                               data->scales.values()[static_cast<std::size_t>(i)])
                            .value;
                physical_phase(data->right.values().subspan(static_cast<std::size_t>(j * n),
                                                            static_cast<std::size_t>(n)),
                               data->scales.values(), PairSide::Right, work.values(),
                               candidate.values());
                if (request.left) {
                    const auto target = n - 1 - diagonal;
                    const auto actual =
                        scaled(data->values.values()[at], 1, data->eigenvalue_scale).value;
                    if (!associated(std::conj(actual), leftt(target, target), request.association))
                        reject(FailureCode::PrecisionUnresolved);
                    eigenvector(leftq, leftt, target, work.values(), candidate.values());
                    for (core::Index i = 0; i < n; ++i)
                        data->left.values()[static_cast<std::size_t>(i + j * n)] =
                            scaled(candidate.values()[static_cast<std::size_t>(i)], 1,
                                   data->scales.values()[static_cast<std::size_t>(i)])
                                .value;
                    physical_phase(data->left.values().subspan(static_cast<std::size_t>(j * n),
                                                               static_cast<std::size_t>(n)),
                                   data->scales.values(), PairSide::Left, work.values(),
                                   candidate.values());
                }
            }
        }
        {
            const FullStageClock qualification{*data, &FullTimings::qualification};
            data->diagnostics.converged_count = static_cast<std::uint64_t>(k);
            if (request.require_full_basis)
                independence(*data, z, inverse);
            qualify_pairs(*data, request, actions.values(), errors.values());
        }
        check_timing(*data);
        return FullPairs{std::move(data)};
    });
}
} // namespace detail
Result<FullPairs> full_pairs(const FullSpectrum& packet, const FullSelection& selected,
                             const FullPairRequest& request, core::ResourcePlan& plan) noexcept {
    return detail::FullFactory::pairs(packet, selected, request, plan);
}
} // namespace uvlas::numerics::eigen
