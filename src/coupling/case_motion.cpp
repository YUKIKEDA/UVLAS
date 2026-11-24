#include "case_motion.hpp"

#include "case_component_precision_private.hpp"

#include <limits>
#include <type_traits>
#include <utility>

namespace uvlas::coupling::detail {
namespace a = aero::detail;
namespace p = a::pchip;
namespace m = motion::m;

namespace {
using enum core::FailureCode;
constexpr auto preparation = core::Operation::PrepareModel;
constexpr auto evaluation = core::Operation::EvaluateStep;
constexpr auto axes = std::size_t{3};

template <class T>
T take(core::Result<T> result) {
    if (!result)
        throw result.error();
    if constexpr (!std::is_void_v<T>)
        return std::move(*result);
}

template <class T, class F>
core::Result<T> guarded(core::Operation operation, core::FailureDetails details,
                        F evaluate) noexcept {
    return core::guard_entry(operation, PrecisionUnresolved, details, [&]() -> core::Result<T> {
        try {
            return evaluate();
        } catch (const core::Failure& failure) {
            return std::unexpected(failure);
        } catch (const m::ArithmeticFailure& failure) {
            return std::unexpected(core::Failure{failure.code, operation, details});
        }
    });
}

constexpr std::array jet_components{&motion::ScalarJet::value, &motion::ScalarJet::first,
                                    &motion::ScalarJet::second};

} // namespace

namespace case_motion_detail {
struct Builder {
    static CaseMotionReady prepare(const io::SiCase& source, CaseMotionField field,
                                   CaseMotionTimeWindow window, p::PchipInputStamp stamp,
                                   p::CurveWorkspace& workspace, const CaseMotionLimits& limits,
                                   core::ResourcePlan& memory, const std::stop_token& stop,
                                   core::FailureDetails details) {
        m::Work work{0, limits.maximum_interval_work, stop};
        work.step();
        auto frame =
            take(CaseMotionFrame::prepare(source, field, limits.frame, memory, stop, details));
        auto curves = take(prepare_case_motion_time_series(
            source, field, window, stamp, workspace, limits.time_series, memory, stop, details));
        auto pi = m::pi_constant(work);
        work.step();
        return CaseMotionReady{std::move(frame), std::move(curves), std::move(pi)};
    }

    static CaseMotionEvaluation evaluate(const CaseMotionReady& prepared,
                                         const CaseMotionReady& current, CaseMotionTimeQuery query,
                                         const CaseMotionProjection& request,
                                         p::CurveWorkspace& workspace,
                                         const CaseMotionLimits& limits, core::ResourcePlan& memory,
                                         const std::stop_token& stop,
                                         core::FailureDetails details) {
        m::Work work{0, limits.maximum_interval_work, stop};
        work.step();
        const auto orders = static_cast<std::size_t>(query.order) + 1;
        constexpr auto width = CaseMotionEvaluation::components_per_order;
        if (orders > jet_components.size() || request.components.size() != orders * width)
            case_components::reject(InvalidInput);
        if (!a::p1_environment_supported() ||
            (request.representation != core::Precision::Binary64 &&
             request.representation != core::Precision::TwoWord))
            case_components::reject(PrecisionUnsupported);
        const auto point_stride = orders * axes;
        constexpr auto max_index =
            static_cast<std::size_t>(std::numeric_limits<core::Index>::max());
        if (request.reference_points.size() > max_index / point_stride)
            case_components::reject(ResourceExceeded);
        for (const auto& point : request.reference_points) {
            work.step();
            if (point.components.size() != point_stride)
                case_components::reject(InvalidInput);
        }
        take(prepared.frame_.match(current.frame_, stop, details));
        const auto frame = take(prepared.frame_.view());
        auto curves = take(evaluate_case_motion_time_series(prepared.curves_, current.curves_,
                                                            query, workspace, limits.time_series,
                                                            memory, stop, details));
        const auto pi_power = curves.rotation_pi_power();
        if (pi_power < 0 || pi_power > 1)
            case_components::reject(InvalidInput);
        motion::VectorJet offset{}, phi{};
        for (std::size_t order = 0; order < orders; ++order) {
            const auto derivative = static_cast<p::JetOrder>(order);
            for (std::size_t axis = 0; axis < axes; ++axis) {
                offset[axis].*jet_components[order] =
                    case_components::rational(take(curves.translation(axis, derivative)),
                                              limits.maximum_interval_input_bits, work);
                auto angle =
                    case_components::rational(take(curves.rotation_coefficient(axis, derivative)),
                                              limits.maximum_interval_input_bits, work);
                phi[axis].*jet_components[order] =
                    pi_power == 0 ? angle : m::multiply(angle, prepared.pi_, work);
            }
        }
        const auto rigid =
            motion::compose_local(frame.origin, frame.rotation, offset, phi, prepared.pi_, work);
        a::LocalRootScratch scratch{};
        CaseMotionEvaluation::Values values{};
        std::size_t slot = 0;
        for (std::size_t order = 0; order < orders; ++order) {
            const auto put = [&](const motion::ScalarJet& jet) {
                values[slot] = case_components::project(
                    jet.*jet_components[order], request.components[slot], request.representation,
                    limits.maximum_interval_input_bits, scratch, work);
                ++slot;
            };
            for (const auto& component : rigid.centre)
                put(component);
            for (const auto& component : rigid.rotation)
                put(component);
        }
        const auto point_words =
            static_cast<core::Index>(request.reference_points.size() * point_stride);
        auto points = take(memory.allocate<io::SiEnclosure>(point_words, evaluation));
        auto output = points.values();
        for (const auto& point : request.reference_points) {
            motion::Vector reference{};
            for (std::size_t axis = 0; axis < axes; ++axis)
                reference[axis] = case_components::rational(
                    point.world_metres[axis], limits.maximum_interval_input_bits, work);
            const auto mapped =
                motion::map_reference_point(reference, frame.origin, frame.rotation, rigid, work);
            std::size_t index = 0;
            for (std::size_t order = 0; order < orders; ++order)
                for (const auto& component : mapped) {
                    output[index] = case_components::project(
                        component.*jet_components[order], point.components[index],
                        request.representation, limits.maximum_interval_input_bits, scratch, work);
                    ++index;
                }
            output = output.subspan(point_stride);
        }
        work.step();
        return CaseMotionEvaluation{std::move(curves), values, std::move(points), orders,
                                    work.used};
    }
};
} // namespace case_motion_detail

CaseMotionReady::CaseMotionReady(CaseMotionFrame frame, CaseMotionTimeSeries curves,
                                 m::Interval pi) noexcept
    : frame_(std::move(frame)), curves_(std::move(curves)), pi_(std::move(pi)) {}

std::span<const p::RationalView> CaseMotionReady::events() const& noexcept {
    return curves_.events();
}

CaseMotionEvaluation::CaseMotionEvaluation(CaseMotionTimeJet time, const Values& values,
                                           core::PlannedBuffer<io::SiEnclosure> points,
                                           std::size_t orders, std::uint64_t work) noexcept
    : time_(std::move(time)),
      values_(values),
      points_(std::move(points)),
      orders_(orders),
      interval_work_(work) {}

CaseMotionEvaluation::CaseMotionEvaluation(CaseMotionEvaluation&& source) noexcept
    : CaseMotionEvaluation(std::move(source.time_), source.values_, std::move(source.points_),
                           std::exchange(source.orders_, 0), source.interval_work_) {}

CaseMotionEvaluation& CaseMotionEvaluation::operator=(CaseMotionEvaluation&& source) noexcept {
    if (this != &source) {
        time_ = std::move(source.time_);
        values_ = source.values_;
        points_ = std::move(source.points_);
        orders_ = std::exchange(source.orders_, 0);
        interval_work_ = source.interval_work_;
    }
    return *this;
}

core::Result<io::SiEnclosure> CaseMotionEvaluation::centre(std::size_t axis,
                                                           p::JetOrder order) const noexcept {
    const auto index = static_cast<std::size_t>(order);
    if (axis >= axes || index >= orders_)
        return std::unexpected(core::Failure{InvalidInput, evaluation});
    return values_[index * components_per_order + axis];
}

core::Result<io::SiEnclosure> CaseMotionEvaluation::rotation(std::size_t row, std::size_t column,
                                                             p::JetOrder order) const noexcept {
    const auto index = static_cast<std::size_t>(order);
    if (row >= axes || column >= axes || index >= orders_)
        return std::unexpected(core::Failure{InvalidInput, evaluation});
    return values_[index * components_per_order + axes + axes * row + column];
}

core::Result<std::reference_wrapper<const CaseMotionTimeJet>>
CaseMotionEvaluation::time() const& noexcept {
    if (orders_ == 0)
        return std::unexpected(core::Failure{InvalidInput, evaluation});
    return std::cref(time_);
}

std::size_t CaseMotionEvaluation::point_count() const noexcept {
    return orders_ == 0 ? 0 : static_cast<std::size_t>(points_.size()) / (axes * orders_);
}

core::Result<io::SiEnclosure> CaseMotionEvaluation::point(std::size_t point, std::size_t axis,
                                                          p::JetOrder order) const noexcept {
    const auto index = static_cast<std::size_t>(order);
    if (point >= point_count() || axis >= axes || index >= orders_)
        return std::unexpected(core::Failure{InvalidInput, evaluation});
    return points_.values()[(point * orders_ + index) * axes + axis];
}

std::uint64_t CaseMotionEvaluation::interval_work() const noexcept {
    return orders_ == 0 ? 0 : interval_work_;
}

core::Result<CaseMotionReady>
prepare_case_motion(const io::SiCase& source, CaseMotionField field, CaseMotionTimeWindow window,
                    p::PchipInputStamp stamp, p::CurveWorkspace& workspace,
                    const CaseMotionLimits& limits, core::ResourcePlan& memory,
                    const std::stop_token& stop, core::FailureDetails details) noexcept {
    return guarded<CaseMotionReady>(preparation, details, [&] {
        return case_motion_detail::Builder::prepare(source, field, window, stamp, workspace, limits,
                                                    memory, stop, details);
    });
}

core::Result<CaseMotionEvaluation>
evaluate_case_motion(const CaseMotionReady& prepared, const CaseMotionReady& current,
                     CaseMotionTimeQuery query, const CaseMotionProjection& request,
                     p::CurveWorkspace& workspace, const CaseMotionLimits& limits,
                     core::ResourcePlan& memory, const std::stop_token& stop,
                     core::FailureDetails details) noexcept {
    return guarded<CaseMotionEvaluation>(evaluation, details, [&] {
        return case_motion_detail::Builder::evaluate(prepared, current, query, request, workspace,
                                                     limits, memory, stop, details);
    });
}

} // namespace uvlas::coupling::detail
