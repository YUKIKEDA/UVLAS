#include "linear_solve_detail.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <utility>

namespace uvlas::numerics::linear {
Certificate::Certificate(
    Context context, core::RhsGeneration rhs, core::NumericalGeneration rhs_buffer,
    Accuracy accuracy, core::ErrorCertificate retained_residual,
    core::ErrorCertificate retained_backward, core::ErrorCertificate residual,
    core::ErrorCertificate backward, core::PlannedBuffer<ColumnMetrics> retained_columns,
    core::PlannedBuffer<ColumnMetrics> columns, core::PlannedBuffer<double> unknown_scales,
    core::PlannedBuffer<double> residual_scales, core::PlannedBuffer<UnitDimensions> unknown_units,
    core::PlannedBuffer<UnitDimensions> residual_units,
    core::PlannedBuffer<std::optional<ConsumerProof>> consumers,
    BackendAdmission preparation_admission, BackendAdmission solve_admission, unsigned refinements,
    RetentionPath retention, ConditionDiagnostic condition,
    LinearSystemOrientation orientation) noexcept
    : context_(context),
      orientation_(orientation),
      rhs_(rhs),
      rhs_buffer_(rhs_buffer),
      accuracy_(accuracy),
      retained_residual_(retained_residual),
      retained_backward_(retained_backward),
      residual_(residual),
      backward_(backward),
      retained_columns_(std::move(retained_columns)),
      columns_(std::move(columns)),
      unknown_scales_(std::move(unknown_scales)),
      residual_scales_(std::move(residual_scales)),
      unknown_units_(std::move(unknown_units)),
      residual_units_(std::move(residual_units)),
      consumers_(std::move(consumers)),
      prepared_admission_(preparation_admission),
      solve_admission_(solve_admission),
      refinements_(refinements),
      retention_(retention),
      condition_(condition) {}

core::Index Certificate::consumer_count() const noexcept {
    // prepare_metadata rejects a count that does not fit the public index domain.
    return static_cast<core::Index>(consumers_.size());
}
core::Result<ConsumerProof> Certificate::consumer_proof(core::Index index) const noexcept {
    if (index < 0 || std::cmp_greater_equal(index, consumers_.size())) {
        return std::unexpected(
            core::Failure{core::FailureCode::InvalidInput, core::Operation::SolveLinear});
    }
    const auto& proof = consumers_.values()[static_cast<std::size_t>(index)];
    if (!proof) {
        return std::unexpected(
            core::Failure{core::FailureCode::PrecisionUnresolved, core::Operation::SolveLinear});
    }
    return *proof;
}

namespace detail {
namespace {
constexpr auto operation = core::Operation::SolveLinear;
auto fail(core::FailureCode code) noexcept {
    return std::unexpected(core::Failure{code, operation});
}
// A failed necessary observation identifies the actual evaluated candidate. Keep any
// more specific diagnostics supplied by the consumer instead of relabelling them.
core::Failure consumer_failure(core::Failure failure, const ConsumerCheck& check,
                               core::NumericalGeneration generation,
                               const Context& context) noexcept {
    auto& details = failure.details;
    if (!details.model) {
        details.model = context.cache.model;
    }
    if (!details.trial) {
        details.trial = context.trial;
    }
    if (!details.numerical_generation) {
        details.numerical_generation = generation;
    }
    if (!details.problem) {
        details.problem = check.observation;
    }
    return failure;
}
bool valid_kind(CheckKind kind) noexcept {
    return kind == CheckKind::Observation || kind == CheckKind::Constraint ||
           kind == CheckKind::Reaction;
}
} // namespace

core::Result<void> check_consumer_inputs(std::span<const ConsumerCheck> checks) noexcept {
    if (!std::in_range<core::Index>(checks.size())) {
        return fail(core::FailureCode::ResourceExceeded);
    }
    for (const auto& check : checks) {
        if (!std::isfinite(check.budget)) {
            return fail(core::FailureCode::Nonfinite);
        }
        if (check.budget < 0 || !valid_kind(check.kind)) {
            return fail(core::FailureCode::InvalidInput);
        }
    }
    return {};
}
core::Result<void> check_accuracy_inputs(const Accuracy& accuracy) noexcept {
    for (const double value : {accuracy.absolute_tolerance, accuracy.relative_tolerance,
                               accuracy.componentwise_backward_tolerance}) {
        if (!std::isfinite(value)) {
            return fail(core::FailureCode::Nonfinite);
        }
        if (value < 0) {
            return fail(core::FailureCode::InvalidInput);
        }
    }
    if (accuracy.norm != Norm::Infinity && accuracy.norm != Norm::L2) {
        return fail(core::FailureCode::InvalidInput);
    }
    return {};
}
namespace {
core::Result<double> aggregate(std::span<const ColumnMetrics> columns, const Accuracy& accuracy,
                               bool backward) noexcept {
    double maximum = 0;
    for (const auto& column : columns) {
        if (!std::isfinite(column.normalized_original_residual_upper) ||
            !std::isfinite(column.componentwise_backward_error_upper) ||
            !std::isfinite(column.scaled_correction_upper)) {
            return fail(core::FailureCode::Nonfinite);
        }
        if (column.normalized_original_residual_upper < 0 ||
            column.componentwise_backward_error_upper < 0 || column.scaled_correction_upper < 0) {
            return fail(core::FailureCode::InvalidInput);
        }
        if (!meets_accuracy(column, accuracy)) {
            return fail(core::FailureCode::LinearUnconverged);
        }
        maximum = std::max(maximum, backward ? column.componentwise_backward_error_upper
                                             : column.normalized_original_residual_upper);
    }
    return maximum;
}
} // namespace

core::Result<MetadataBuffers>
SolutionFactory::prepare_metadata(core::Index rows, core::Index rhs_columns, const Scales& scales,
                                  std::span<const ConsumerCheck> checks,
                                  core::ResourcePlan& result_plan) noexcept {
    if (rows < 0 || rhs_columns < 0 || std::cmp_not_equal(scales.unknown.size(), rows) ||
        std::cmp_not_equal(scales.residual.size(), rows) ||
        std::cmp_not_equal(scales.unknown_units.size(), rows) ||
        std::cmp_not_equal(scales.residual_units.size(), rows)) {
        return fail(core::FailureCode::InvalidInput);
    }
    if (!std::isfinite(scales.context.time)) {
        return fail(core::FailureCode::Nonfinite);
    }
    for (std::size_t i = 0; i < scales.unknown.size(); ++i) {
        if (!std::isfinite(scales.unknown[i]) || !std::isfinite(scales.residual[i])) {
            return fail(core::FailureCode::Nonfinite);
        }
        if (scales.unknown[i] <= 0 || scales.residual[i] <= 0) {
            return fail(core::FailureCode::InvalidInput);
        }
    }
    const auto valid = check_consumer_inputs(checks);
    if (!valid) {
        return std::unexpected(valid.error());
    }
    auto retained = result_plan.allocate<ColumnMetrics>(rhs_columns, operation);
    if (!retained) {
        return std::unexpected(retained.error());
    }
    auto returned = result_plan.allocate<ColumnMetrics>(rhs_columns, operation);
    if (!returned) {
        return std::unexpected(returned.error());
    }
    auto unknown = result_plan.allocate<double>(rows, operation);
    if (!unknown) {
        return std::unexpected(unknown.error());
    }
    auto residual = result_plan.allocate<double>(rows, operation);
    if (!residual) {
        return std::unexpected(residual.error());
    }
    auto unknown_units = result_plan.allocate<UnitDimensions>(rows, operation);
    if (!unknown_units) {
        return std::unexpected(unknown_units.error());
    }
    auto residual_units = result_plan.allocate<UnitDimensions>(rows, operation);
    if (!residual_units) {
        return std::unexpected(residual_units.error());
    }
    auto consumers = result_plan.allocate<std::optional<ConsumerProof>>(
        static_cast<core::Index>(checks.size()), operation);
    if (!consumers) {
        return std::unexpected(consumers.error());
    }
    // An unwritten column is unevaluated, never an exact-zero certificate.
    const double unset = std::numeric_limits<double>::quiet_NaN();
    std::ranges::fill(retained->values(), ColumnMetrics{unset, unset, unset});
    std::ranges::fill(returned->values(), ColumnMetrics{unset, unset, unset});
    std::ranges::copy(scales.unknown, unknown->values().begin());
    std::ranges::copy(scales.residual, residual->values().begin());
    std::ranges::copy(scales.unknown_units, unknown_units->values().begin());
    std::ranges::copy(scales.residual_units, residual_units->values().begin());
    return MetadataBuffers{std::move(*retained),      std::move(*returned),
                           std::move(*unknown),       std::move(*residual),
                           std::move(*unknown_units), std::move(*residual_units),
                           std::move(*consumers)};
}

template <StorageScalar Scalar>
core::Result<Solution<Scalar>> SolutionFactory::finalize(
    PrecisionBuffer<Scalar> retained, std::optional<PrecisionBuffer<Scalar>> projected,
    const Precision<Scalar>& precision, MetadataBuffers metadata, Finalization context) noexcept {
    const auto valid_accuracy = check_accuracy_inputs(context.accuracy);
    if (!valid_accuracy) {
        return std::unexpected(valid_accuracy.error());
    }
    const auto valid_checks = check_consumer_inputs(precision.required_checks);
    if (!valid_checks) {
        return std::unexpected(valid_checks.error());
    }
    const auto shape = retained.shape();
    if (shape.rows < 0 || shape.columns < 0 ||
        std::cmp_not_equal(metadata.retained_columns.size(), shape.columns) ||
        std::cmp_not_equal(metadata.returned_columns.size(), shape.columns) ||
        std::cmp_not_equal(metadata.unknown_scales.size(), shape.rows) ||
        std::cmp_not_equal(metadata.residual_scales.size(), shape.rows) ||
        std::cmp_not_equal(metadata.unknown_units.size(), shape.rows) ||
        std::cmp_not_equal(metadata.residual_units.size(), shape.rows) ||
        metadata.consumers.size() != precision.required_checks.size() || context.refinements > 5) {
        return fail(core::FailureCode::InvalidInput);
    }
    if (!std::isfinite(context.context.time)) {
        return fail(core::FailureCode::Nonfinite);
    }
    const bool two_word = retained.precision() == core::Precision::TwoWord;
    if ((context.retention == RetentionPath::Binary64 && two_word) ||
        ((context.retention == RetentionPath::InitiallyTwoWord ||
          context.retention == RetentionPath::PromotedTwoWord) &&
         !two_word) ||
        (context.retention != RetentionPath::Binary64 &&
         context.retention != RetentionPath::InitiallyTwoWord &&
         context.retention != RetentionPath::PromotedTwoWord)) {
        return fail(core::FailureCode::InvalidInput);
    }
    if (const auto* diagnostic = std::get_if<ReciprocalOneNormEstimate>(&context.condition)) {
        if (!std::isfinite(diagnostic->value)) {
            return fail(core::FailureCode::Nonfinite);
        }
        if (diagnostic->value < 0) {
            return fail(core::FailureCode::InvalidInput);
        }
    }
    if (const auto* request = std::get_if<RetainedPrecision>(&precision.returned)) {
        if (projected) {
            return fail(core::FailureCode::InvalidInput);
        }
        if (request->generation != retained.generation()) {
            return fail(core::FailureCode::StaleGeneration);
        }
        std::ranges::copy(metadata.retained_columns.values(),
                          metadata.returned_columns.values().begin());
    } else {
        const auto* projection_request = std::get_if<Binary64Projection>(&precision.returned);
        if (projection_request == nullptr) {
            return fail(core::FailureCode::PrecisionUnsupported);
        }
        if (!projected || projected->precision() != core::Precision::Binary64) {
            return fail(core::FailureCode::PrecisionUnsupported);
        }
        const auto target_shape = projected->shape();
        if (target_shape.rows != shape.rows || target_shape.columns != shape.columns ||
            target_shape.leading_dimension != shape.leading_dimension) {
            return fail(core::FailureCode::InvalidInput);
        }
        if (retained.generation() != projection_request->retained_generation ||
            projected->generation() != projection_request->target_generation ||
            projection_request->retained_generation == projection_request->target_generation) {
            return fail(core::FailureCode::StaleGeneration);
        }
    }
    const auto rr = aggregate(metadata.retained_columns.values(), context.accuracy, false);
    if (!rr) {
        return std::unexpected(rr.error());
    }
    const auto rb = aggregate(metadata.retained_columns.values(), context.accuracy, true);
    if (!rb) {
        return std::unexpected(rb.error());
    }
    const auto cr = aggregate(metadata.returned_columns.values(), context.accuracy, false);
    if (!cr) {
        return std::unexpected(cr.error());
    }
    const auto cb = aggregate(metadata.returned_columns.values(), context.accuracy, true);
    if (!cb) {
        return std::unexpected(cb.error());
    }
    const auto& candidate = projected ? *projected : retained;
    const core::CertificateContext retained_residual_context{
        context.accuracy.residual_observation, context.context.residual_layout,
        retained.generation(), retained.precision()};
    const core::CertificateContext retained_backward_context{
        context.accuracy.backward_observation, context.context.residual_layout,
        retained.generation(), retained.precision()};
    const core::CertificateContext returned_residual_context{
        context.accuracy.residual_observation, context.context.residual_layout,
        candidate.generation(), candidate.precision()};
    const core::CertificateContext returned_backward_context{
        context.accuracy.backward_observation, context.context.residual_layout,
        candidate.generation(), candidate.precision()};
    auto retained_residual = core::ErrorCertificate::assess(retained_residual_context, 1, *rr);
    auto retained_backward = core::ErrorCertificate::assess(
        retained_backward_context, context.accuracy.componentwise_backward_tolerance, *rb);
    auto residual = core::ErrorCertificate::assess(returned_residual_context, 1, *cr);
    auto backward = core::ErrorCertificate::assess(
        returned_backward_context, context.accuracy.componentwise_backward_tolerance, *cb);
    if (!retained_residual) {
        return std::unexpected(retained_residual.error());
    }
    if (!retained_backward) {
        return std::unexpected(retained_backward.error());
    }
    if (!residual) {
        return std::unexpected(residual.error());
    }
    if (!backward) {
        return std::unexpected(backward.error());
    }

    // An algebraic certificate never substitutes for a required physical/constraint/reaction check.
    for (std::size_t i = 0; i < precision.required_checks.size(); ++i) {
        const auto& check = precision.required_checks[i];
        auto supported =
            require_precision(retained.precision(), precision.evaluator.supported_precision());
        if (!supported) {
            return std::unexpected(
                consumer_failure(supported.error(), check, retained.generation(), context.context));
        }
        supported =
            require_precision(candidate.precision(), precision.evaluator.supported_precision());
        if (!supported) {
            return std::unexpected(consumer_failure(supported.error(), check,
                                                    candidate.generation(), context.context));
        }
        const core::CertificateContext expected{check.observation, check.scope,
                                                candidate.generation(), candidate.precision()};
        auto proof =
            core::guard_entry(operation, core::FailureCode::PrecisionUnresolved, {}, [&]() {
                return precision.evaluator.assess(retained, candidate, context.context, check,
                                                  expected);
            });
        if (!proof) {
            return std::unexpected(
                consumer_failure(proof.error(), check, candidate.generation(), context.context));
        }
        if (proof->context() != expected || std::bit_cast<std::uint64_t>(proof->budget()) !=
                                                std::bit_cast<std::uint64_t>(check.budget)) {
            return std::unexpected(
                consumer_failure(core::Failure{core::FailureCode::StaleGeneration, operation},
                                 check, candidate.generation(), context.context));
        }
        if (!std::holds_alternative<core::Achieved>(proof->assessment())) {
            return std::unexpected(
                consumer_failure(core::Failure{core::FailureCode::PrecisionUnresolved, operation},
                                 check, candidate.generation(), context.context));
        }
        metadata.consumers.values()[i].emplace(ConsumerProof{check, *proof});
    }
    Certificate certificate{context.context,
                            context.rhs,
                            context.rhs_buffer_generation,
                            context.accuracy,
                            *retained_residual,
                            *retained_backward,
                            *residual,
                            *backward,
                            std::move(metadata.retained_columns),
                            std::move(metadata.returned_columns),
                            std::move(metadata.unknown_scales),
                            std::move(metadata.residual_scales),
                            std::move(metadata.unknown_units),
                            std::move(metadata.residual_units),
                            std::move(metadata.consumers),
                            context.preparation_admission,
                            context.solve_admission,
                            context.refinements,
                            context.retention,
                            context.condition,
                            context.orientation};
    return Solution<Scalar>{projected ? std::move(*projected) : std::move(retained),
                            std::move(certificate)};
}

template core::Result<Solution<double>>
SolutionFactory::finalize(PrecisionBuffer<double>, std::optional<PrecisionBuffer<double>>,
                          const Precision<double>&, MetadataBuffers, Finalization) noexcept;
template core::Result<Solution<std::complex<double>>> SolutionFactory::finalize(
    PrecisionBuffer<std::complex<double>>, std::optional<PrecisionBuffer<std::complex<double>>>,
    const Precision<std::complex<double>>&, MetadataBuffers, Finalization) noexcept;
} // namespace detail
} // namespace uvlas::numerics::linear
