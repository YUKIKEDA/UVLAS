#include "case_uniform_flow_projection.hpp"

#include "io/si_rational_projection.hpp"

namespace uvlas::coupling::detail {
namespace {
constexpr std::size_t components = 3;
constexpr auto operation = core::Operation::EvaluateAero;
constexpr auto invalid = core::FailureCode::InvalidInput;

io::detail::SiRationalWordsView words(aero::detail::pchip::RationalView value) noexcept {
    return {value.sign, value.numerator, value.denominator};
}
io::detail::SiRationalProjectionBudget budget(const CaseFlowComponentBudget& component) noexcept {
    return {words(component.scale), words(component.absolute), words(component.relative)};
}
} // namespace

ProjectedCaseUniformFlow::ProjectedCaseUniformFlow(CaseUniformFlowEvaluation&& parent,
                                                   const std::array<io::SiEnclosure, 6>& values,
                                                   core::Precision representation,
                                                   std::uint64_t work) noexcept
    : source_(std::move(parent)), values_(values), representation_(representation), work_(work) {}

const CaseUniformFlowEvaluation& ProjectedCaseUniformFlow::source() const& noexcept {
    return source_;
}

core::Result<io::SiEnclosure>
ProjectedCaseUniformFlow::velocity(std::size_t component) const noexcept {
    if (!source_.jet.active() || component >= components)
        return std::unexpected(core::Failure{invalid, operation});
    return values_[component];
}
core::Result<io::SiEnclosure>
ProjectedCaseUniformFlow::time_partial(std::size_t component) const noexcept {
    if (!source_.jet.active() || component >= components)
        return std::unexpected(core::Failure{invalid, operation});
    return values_[components + component];
}
core::Result<io::SiEnclosure>
ProjectedCaseUniformFlow::spatial_partial(std::size_t component,
                                          std::size_t coordinate) const noexcept {
    if (!source_.jet.active() || component >= components || coordinate >= components)
        return std::unexpected(core::Failure{invalid, operation});
    if (representation_ == core::Precision::TwoWord)
        return io::SiEnclosure{io::ScaledTwoWord{0.0, 0.0, 0.0, 0}};
    return io::SiEnclosure{io::ScaledBinary64{0.0, 0.0, 0}};
}
core::Precision ProjectedCaseUniformFlow::representation() const noexcept {
    return representation_;
}
std::uint64_t ProjectedCaseUniformFlow::used_steps() const noexcept {
    return source_.jet.active() ? work_ : 0;
}

core::Result<ProjectedCaseUniformFlow>
project_case_uniform_flow(CaseUniformFlowEvaluation&& source,
                          const CaseFlowProjectionRequest& request, const std::stop_token& stop,
                          core::FailureDetails details) noexcept {
    if (!source.jet.active())
        return std::unexpected(core::Failure{invalid, operation, details});
    io::detail::SiRationalProjectionLimits remaining{request.limits.input_integer_bits,
                                                     request.limits.steps};
    std::array<io::SiEnclosure, 2 * components> values{};
    for (std::size_t index = 0; index < values.size(); ++index) {
        const bool velocity = index < components;
        const auto component = index % components;
        const auto exact =
            velocity ? source.jet.velocity(component) : source.jet.time_partial(component);
        if (!exact)
            return std::unexpected(exact.error());
        const auto& accuracy =
            velocity ? request.velocity[component] : request.time_partial[component];
        const auto candidate = io::detail::project_rational_words(
            words(*exact), budget(accuracy), request.representation, remaining, stop, details);
        if (!candidate)
            return std::unexpected(candidate.error());
        values[index] = candidate->value;
        remaining.steps -= candidate->used_steps;
    }
    if (stop.stop_requested())
        return std::unexpected(core::Failure{core::FailureCode::Cancelled, operation, details});
    return ProjectedCaseUniformFlow{std::move(source), values, request.representation,
                                    request.limits.steps - remaining.steps};
}

namespace case_flow_detail {
struct PrescribedProjection {
    static core::Result<ProjectedCasePrescribedWake>
    project(CasePrescribedWakeEvaluation&& source, const CaseFlowProjectionRequest& request,
            const std::stop_token& stop, core::FailureDetails details) noexcept {
        auto result = project_case_uniform_flow(std::move(source.value_), request, stop, details);
        if (!result)
            return std::unexpected(result.error());
        return ProjectedCasePrescribedWake{std::move(*result)};
    }
};
} // namespace case_flow_detail

core::Result<ProjectedCasePrescribedWake>
project_case_prescribed_wake(CasePrescribedWakeEvaluation&& source,
                             const CaseFlowProjectionRequest& request, const std::stop_token& stop,
                             core::FailureDetails details) noexcept {
    return case_flow_detail::PrescribedProjection::project(std::move(source), request, stop,
                                                           details);
}

} // namespace uvlas::coupling::detail
