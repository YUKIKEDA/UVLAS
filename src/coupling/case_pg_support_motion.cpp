#include "case_pg_support_motion.hpp"

#include "case_component_precision_private.hpp"
#include "structure/exact_polynomial_detail.hpp"

#include <limits>

namespace uvlas::coupling::detail {
namespace {
namespace a = aero::detail;
namespace p = a::pchip;
namespace m = a::p1;
using enum core::FailureCode;
constexpr auto operation = core::Operation::EvaluateStep;
constexpr std::size_t axes = 3;

template <class T>
T take(core::Result<T> result) {
    if (!result)
        throw result.error();
    return std::move(*result);
}

m::Interval rational(const structure::ExactRational& value, std::uint64_t limit, m::Work& work) {
    work.step();
    constexpr auto storage = std::numeric_limits<a::ExactInteger>::digits;
    if (value.numerator_bits() > limit || value.denominator_bits() > limit ||
        value.numerator_bits() > storage || value.denominator_bits() > storage)
        case_components::reject(ResourceExceeded);
    if (value.sign() == 0)
        return m::integer(0);

    using n = structure::detail::IntegerBackend;
    const auto assemble = [&](const n& source) {
        a::ExactInteger result = 0;
        for (std::size_t i = source.size(); i > 0;) {
            work.step();
            --i;
            result <<= n::limb_bits;
            result += source.limbs()[i];
        }
        return source.sign() ? -result : result;
    };
    const auto& rep = structure::detail::ExactFactory::value(value).backend();
    return m::divide(m::exact(assemble(rep.num()), 0), m::exact(assemble(rep.denom()), 0), work);
}

m::Interval interval(const io::SiEnclosure& value, m::Work& work) {
    return std::visit(
        [&](const auto& word) {
            m::Interval centre;
            if constexpr (std::is_same_v<std::remove_cvref_t<decltype(word)>, io::ScaledBinary64>)
                centre = a::decode_p1_word(word.center);
            else
                centre = m::add(a::decode_p1_word(word.high), a::decode_p1_word(word.low), work);
            const auto exponent = m::exponent(word.exponent2);
            const auto radius = a::decode_p1_word(word.radius);
            const auto lower = m::scale(m::subtract(centre, radius, work), exponent);
            const auto upper = m::scale(m::add(centre, radius, work), exponent);
            return m::Interval{lower.lower, upper.upper};
        },
        value);
}

std::array<CaseMotionReferencePoint, 4>
reference_probes(std::span<const CaseFlowComponentBudget> budgets) {
    static constexpr std::array<std::uint32_t, 1> one{1};
    const p::RationalView zero{0, {}, one}, basis{1, one, one};
    return {{{{zero, zero, zero}, budgets},
             {{basis, zero, zero}, budgets},
             {{zero, basis, zero}, budgets},
             {{zero, zero, basis}, budgets}}};
}

struct RigidCoefficients {
    std::array<m::Interval, axes> translation;
    std::array<std::array<m::Interval, axes>, axes> rotation;
};

RigidCoefficients coefficients(const CaseMotionEvaluation& motion, p::JetOrder order,
                               m::Work& work) {
    RigidCoefficients out{};
    for (std::size_t axis = 0; axis < axes; ++axis) {
        out.translation[axis] = interval(take(motion.point(0, axis, order)), work);
        for (std::size_t j = 0; j < axes; ++j) {
            const auto mapped = interval(take(motion.point(j + 1, axis, order)), work);
            out.rotation[axis][j] = m::subtract(mapped, out.translation[axis], work);
        }
    }
    return out;
}

void fill_rhs(const structure::PgRigidSupportRows& support, const CaseMotionEvaluation& motion,
              const CasePgSupportProjection& projection, const CasePgSupportLimits& limits,
              std::size_t orders, std::span<io::SiEnclosure> output, m::Work& work) {
    std::array<RigidCoefficients, 3> jets{};
    for (std::size_t d = 0; d < orders; ++d)
        jets[d] = coefficients(motion, static_cast<p::JetOrder>(d), work);
    auto scratch = a::LocalRootScratch{};
    const auto map = support.rows();
    const auto original = support.chart().original().original_rows();
    for (std::size_t row = 0; row < map.size(); ++row) {
        work.step();
        const auto& source = map[row];
        std::array<m::Interval, axes> x{};
        for (std::size_t k = 0; k < axes; ++k)
            x[k] =
                rational(source.moving ? source.reference_position[k] : original[row].prescribed[k],
                         limits.maximum_row_width_bits, work);
        const auto b = source.moving
                           ? rational(source.translation, limits.maximum_row_width_bits, work)
                           : m::integer(0);
        for (std::size_t d = 0; d < orders; ++d)
            for (std::size_t axis = 0; axis < axes; ++axis) {
                auto value = d == 0 ? x[axis] : m::integer(0);
                if (source.moving) {
                    value = m::multiply(jets[d].translation[axis], b, work);
                    for (std::size_t k = 0; k < axes; ++k)
                        value =
                            m::add(value, m::multiply(jets[d].rotation[axis][k], x[k], work), work);
                }
                output[(d * map.size() + row) * axes + axis] = case_components::project(
                    value, projection.original_rhs[d * axes + axis], projection.representation,
                    limits.maximum_row_width_bits, scratch, work);
            }
    }
}
} // namespace

core::Result<CasePgSupportSample> evaluate_case_clocked_pg_support(
    const structure::PgRigidSupportRows& support, const CaseClockedPrescribedInputs& field,
    const CaseClockedPrescribedInputs& current, CaseFlowClockQuery address,
    const std::array<p::RationalView, 3>& flow_position, const CaseKnownFlowProjection& flow,
    const std::optional<CaseFlowProjectionRequest>& wake, p::JetOrder order,
    const CasePgSupportProjection& projection, a::UniformFlowWorkspace& flow_workspace,
    p::CurveWorkspace& motion_workspace, const CasePgSupportLimits& limits,
    core::ResourcePlan& plan, const std::stop_token& stop, core::FailureDetails details) noexcept {
    return core::guard_entry(
        operation, PrecisionUnresolved, details, [&]() -> core::Result<CasePgSupportSample> {
            try {
                m::Work work{0, limits.maximum_row_work, stop};
                work.step();
                const auto orders = static_cast<std::size_t>(order) + 1;
                if (orders > 3 || projection.motion.size() != 12 * orders ||
                    projection.probe.size() != axes * orders ||
                    projection.original_rhs.size() != axes * orders)
                    case_components::reject(InvalidInput);
                const auto identity = take(support.chart().identity());
                const auto count = support.rows().size();
                if (count > limits.maximum_original_rows ||
                    count > static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max()) /
                                (axes * orders))
                    case_components::reject(ResourceExceeded);
                const auto entries = static_cast<core::Index>(count * axes * orders);
                if (take(core::checked_byte_count(entries, sizeof(io::SiEnclosure), operation)) >
                    plan.remaining_bytes())
                    case_components::reject(ResourceExceeded);
                const auto probes = reference_probes(projection.probe);
                auto sample = take(evaluate_case_clocked_prescribed_inputs(
                    field, current, address, flow_position, flow, wake, order,
                    {projection.representation, projection.motion, probes}, flow_workspace,
                    motion_workspace, limits.clock, plan, stop, details));
                auto rhs = take(plan.allocate<io::SiEnclosure>(entries, operation));
                fill_rhs(support, sample.motion, projection, limits, orders, rhs.values(), work);
                work.step();
                return CasePgSupportSample{
                    std::move(sample), identity, std::move(rhs), count, orders, work.used};
            } catch (const core::Failure& failure) {
                return std::unexpected(failure);
            } catch (const m::ArithmeticFailure& failure) {
                return std::unexpected(core::Failure{failure.code, operation, details});
            }
        });
}

} // namespace uvlas::coupling::detail
