#include "case_aero_attachment.hpp"

#include "aero/query_coordinates_private.hpp"
#include "case_component_precision_private.hpp"
#include "case_dyadic_words_private.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace uvlas::coupling::detail {
namespace a = aero::detail;
namespace p = a::pchip;
namespace m = a::p1;

namespace {
using enum core::FailureCode;
constexpr auto evaluation = core::Operation::EvaluateStep;
constexpr std::size_t axes = 3;

template <class T>
T take(core::Result<T> value) {
    if (!value)
        throw value.error();
    return std::move(*value);
}
template <class T>
T take(a::AeroResult<T> value) {
    if (!value)
        throw value.error().cause;
    return std::move(*value);
}
template <class T>
T take(a::GlobalResult<T> value) {
    if (!value)
        throw value.error().origin.cause;
    return std::move(*value);
}
[[noreturn]] void reject(core::FailureCode code) {
    throw m::ArithmeticFailure{code};
}

void charge(m::Work& work, std::uint64_t quantity) {
    if (quantity > work.limit - work.used)
        reject(ResourceExceeded);
    work.used += quantity;
}
template <class T>
core::PlannedBuffer<T> allocate(core::ResourcePlan& memory, std::size_t count) {
    if (count > static_cast<std::size_t>(std::numeric_limits<core::Index>::max()))
        reject(ResourceExceeded);
    return take(memory.allocate<T>(static_cast<core::Index>(count), evaluation));
}

m::Interval interval(a::ScalarWords words, m::Work& work) {
    if (const auto* one = std::get_if<a::SingleWord>(&words))
        return a::decode_p1_word(one->value);
    const auto& both = std::get<a::DoubleWord>(words);
    return m::add(a::decode_p1_word(both.high), a::decode_p1_word(both.low), work);
}

double unscale(double word, std::int64_t exponent) {
    if (word == 0)
        return word;
    if (exponent < -2200 || exponent > 2200)
        reject(PrecisionUnresolved);
    const auto shift = static_cast<int>(exponent);
    const auto value = std::scalbn(word, shift);
    if (!std::isfinite(value) || std::scalbn(value, -shift) != word)
        reject(PrecisionUnresolved);
    return value;
}
a::ScalarWords centre(io::SiEnclosure value) {
    if (const auto* one = std::get_if<io::ScaledBinary64>(&value))
        return a::SingleWord{unscale(one->center, one->exponent2)};
    const auto& both = std::get<io::ScaledTwoWord>(value);
    return a::DoubleWord{unscale(both.high, both.exponent2), unscale(both.low, both.exponent2)};
}

void qualify_point(m::Interval value, const io::SiEnclosure& reference,
                   const CaseFlowComponentBudget& budget, std::uint64_t bits, m::Work& work) {
    const auto representative = interval(centre(reference), work);
    const auto envelope = std::visit(
        [&](const auto& word) {
            const auto radius =
                m::scale(a::decode_p1_word(word.radius), m::exponent(word.exponent2));
            return m::add(representative,
                          {{-radius.upper.mantissa, radius.upper.power}, radius.upper}, work);
        },
        reference);
    const auto scale = case_components::rational(budget.scale, bits, work);
    const auto absolute = case_components::rational(budget.absolute, bits, work);
    const auto relative = case_components::rational(budget.relative, bits, work);
    const auto tolerance = m::add(m::multiply(scale, absolute, work),
                                  m::multiply(relative, m::absolute(envelope), work), work);
    const auto error = m::absolute(m::subtract(value, envelope, work));
    if (m::compare(error.upper, tolerance.lower) > 0)
        reject(PrecisionUnresolved);
}

std::size_t moving_index(std::span<const a::PointId> moving, a::PointId owner, m::Work& work) {
    const auto at = std::ranges::lower_bound(moving, owner, [&](auto lhs, auto rhs) {
        work.step();
        return lhs < rhs;
    });
    if (at == moving.end() || *at != owner)
        return moving.size();
    return static_cast<std::size_t>(at - moving.begin());
}

std::array<NonpenetrationVector, 2> tangents(const a::OwnedGeometry& geometry, a::FacetId facet,
                                             a::LocalRootScratch& scratch, m::Work& work) {
    const auto ids = take(geometry.facet_points(facet));
    const auto origin = take(geometry.point_words(ids[0], a::Configuration::Current));
    std::array<NonpenetrationVector, 2> result{
        NonpenetrationVector{a::SingleWord{0}, a::SingleWord{0}, a::SingleWord{0}},
        NonpenetrationVector{a::SingleWord{0}, a::SingleWord{0}, a::SingleWord{0}}};
    for (std::size_t t = 0; t < result.size(); ++t) {
        const auto endpoint = take(geometry.point_words(ids[t + 1], a::Configuration::Current));
        for (std::size_t k = 0; k < axes; ++k) {
            const auto diff =
                m::subtract(interval(endpoint[k], work), interval(origin[k], work), work);
            const auto precise = a::enclose_p1_value(diff, core::Precision::TwoWord, scratch);
            const auto words = precise.center();
            const auto* pair = std::get_if<a::DoubleWord>(&words);
            if (precise.radius() != 0 || pair == nullptr)
                reject(PrecisionUnresolved);
            result[t][k] = a::DoubleWord{unscale(pair->high, precise.exponent()),
                                         unscale(pair->low, precise.exponent())};
        }
    }
    return result;
}
} // namespace

namespace aero_attachment_detail {
struct Row {
    a::ResolvedQuery point;
    std::array<NonpenetrationVector, 2> surface;
    NonpenetrationVector velocity;
    PrescribedAeroRow source;
};
struct PointOrdinal {
    a::PointId id{0};
    std::size_t column{};
};
struct Trial {
    std::optional<CaseMotionEvaluation> motion;
    core::PlannedBuffer<PointOrdinal> points;
    std::size_t moving_count{};
    std::optional<a::GlobalContactState> contact;
    core::PlannedBuffer<std::optional<Row>> rows;
};

struct Builder {
    static PrescribedAeroAttachment
    prepare(const io::SiCase& source, const io::schema::AeroAttachmentPrescribed& attachment,
            CaseMotionTimeWindow window, p::PchipInputStamp stamp, p::CurveWorkspace& curves,
            const CaseMotionLimits& limits, core::ResourcePlan& memory, const std::stop_token& stop,
            core::FailureDetails details) {
        return PrescribedAeroAttachment{take(prepare_case_motion(
            source, attachment.motion, window, stamp, curves, limits, memory, stop, details))};
    }

    static PrescribedAeroTrial
    evaluate(const PrescribedAeroAttachment& prepared, const PrescribedAeroAttachment& current,
             const PrescribedAeroReference& reference, CaseMotionTimeQuery query,
             const PrescribedAeroProjection& projection, p::CurveWorkspace& curves,
             a::GeometryWorkspace& coordinates, a::GlobalWorkspace& global,
             const PrescribedAeroLimits& limits, core::ResourcePlan& memory,
             const std::stop_token& stop, core::FailureDetails details) {
        m::Work work{0, limits.work, stop};
        work.step();
        const auto& moving = reference.moving_points;
        const auto& original = reference.geometry;
        const auto& rows = reference.rows;
        if (query.order != p::JetOrder::First || moving.empty() || rows.empty() ||
            projection.point.size() != 2 * axes)
            reject(InvalidInput);
        const auto source = take(reference.topology.geometry());
        const auto& base = source.get();
        constexpr auto largest = static_cast<std::size_t>(std::numeric_limits<core::Index>::max());
        if (original.points.size() > largest / axes || rows.size() > largest ||
            moving.size() > largest - rows.size())
            reject(ResourceExceeded);
        auto packed = allocate<std::array<DyadicWords, axes>>(memory, moving.size() + rows.size());
        auto requests = allocate<CaseMotionReferencePoint>(memory, packed.size());
        const auto resolve = [&](const a::OwnedGeometry& owner, const a::QueryRequest& request) {
            auto budget = limits.geometry;
            budget.max_steps = std::min(budget.max_steps, work.limit - work.used);
            auto result = take(a::resolve_query(owner, request, owner.context(), budget,
                                                coordinates, memory, stop));
            charge(work, budget.max_steps);
            return result;
        };
        const auto position = [&](const a::OwnedGeometry& owner, const a::ResolvedQuery& point,
                                  a::Configuration config) {
            auto budget = limits.geometry;
            budget.max_steps = std::min(budget.max_steps, work.limit - work.used);
            auto result =
                take(a::exact_query_coordinates(owner, point, budget, coordinates, stop, config));
            charge(work, result.used_steps);
            return result;
        };
        const auto prepare_point = [&](std::size_t i, const a::QueryRequest& point) {
            const auto resolved = resolve(base, point);
            const auto values = position(base, resolved, a::Configuration::Reference);
            auto& xyz = packed.values()[i];
            for (std::size_t k = 0; k < axes; ++k)
                xyz[k] = DyadicWords{values.numerator[k],
                                     a::ExactQueryCoordinates::denominator_shift, work};
            requests.values()[i] = {{xyz[0].view(), xyz[1].view(), xyz[2].view()},
                                    projection.point};
        };
        for (std::size_t i = 0; i < moving.size(); ++i) {
            work.step();
            if ((i != 0 && moving[i - 1] >= moving[i]) || take(base.owner(moving[i])) != moving[i])
                reject(InvalidInput);
            prepare_point(i, {a::QueryId{moving[i].value()}, a::SourceVertex{moving[i]}});
        }
        for (std::size_t i = 0; i < rows.size(); ++i) {
            work.step();
            const auto& row = rows[i];
            for (const auto vertex : take(base.facet_points(row.material_point.facet))) {
                if (moving_index(moving, take(base.owner(vertex)), work) == moving.size())
                    reject(InvalidInput);
            }
            prepare_point(moving.size() + i, {row.query, row.material_point});
        }

        auto motion = take(
            evaluate_case_motion(prepared.motion_, current.motion_, query,
                                 {projection.representation, projection.pose, requests.values()},
                                 curves, limits.motion, memory, stop, details));
        charge(work, motion.interval_work());
        const auto instant = take(motion.time()).get().time();
        if (const auto& expected_time = original.exact_time) {
            const auto equal = [&](std::uint32_t left, std::uint32_t right) {
                work.step();
                return left == right;
            };
            if (expected_time->sign != instant.sign ||
                !std::ranges::equal(expected_time->numerator, instant.numerator, equal) ||
                !std::ranges::equal(expected_time->denominator, instant.denominator, equal))
                reject(StaleGeneration);
        }
        const auto word_count = axes * original.points.size();
        auto high = allocate<double>(memory, word_count);
        core::PlannedBuffer<double> low;
        const auto* one = std::get_if<a::SinglePositions>(&original.current);
        if (projection.representation == core::Precision::Binary64 && one != nullptr) {
            if (one->words.size() != word_count)
                reject(InvalidInput);
            std::ranges::copy(one->words, high.values().begin());
        } else if (projection.representation == core::Precision::TwoWord && one == nullptr) {
            const auto& both = std::get<a::DoublePositions>(original.current);
            if (both.high.size() != word_count || both.low.size() != word_count)
                reject(InvalidInput);
            low = allocate<double>(memory, word_count);
            std::ranges::copy(both.high, high.values().begin());
            std::ranges::copy(both.low, low.values().begin());
        } else {
            reject(PrecisionUnsupported);
        }
        charge(work, word_count);
        auto point_rows = allocate<PointOrdinal>(memory, original.points.size());
        for (std::size_t i = 0; i < original.points.size(); ++i) {
            const auto owner = take(base.owner(original.points[i]));
            const auto column = moving_index(moving, owner, work);
            point_rows.values()[i] = {original.points[i], column};
            if (column == moving.size())
                continue;
            for (std::size_t k = 0; k < axes; ++k) {
                work.step();
                const auto words = centre(take(motion.point(column, k, p::JetOrder::Value)));
                if (const auto* word = std::get_if<a::SingleWord>(&words))
                    high.values()[axes * i + k] = word->value;
                else {
                    const auto& pair = std::get<a::DoubleWord>(words);
                    high.values()[axes * i + k] = pair.high;
                    low.values()[axes * i + k] = pair.low;
                }
            }
        }
        auto candidate = original;
        if (low.size() == 0)
            candidate.current = a::SinglePositions{high.values()};
        else
            candidate.current = a::DoublePositions{high.values(), low.values()};
        candidate.exact_time = instant;
        auto geometry_budget = limits.geometry;
        geometry_budget.max_steps = std::min(geometry_budget.max_steps, work.limit - work.used);
        const auto geometry =
            take(a::prepare_owned_geometry(candidate, geometry_budget, coordinates, memory, stop));
        charge(work, geometry_budget.max_steps);
        auto contact_budget = limits.contact;
        contact_budget.max_steps = std::min(contact_budget.max_steps, work.limit - work.used);
        auto contact = take(a::check_global_geometry(reference.topology, geometry, contact_budget,
                                                     global, memory, stop));
        charge(work, take(contact.work_plan()).executed_steps);
        auto owner = allocate<Trial>(memory, 1);
        auto& out = owner.values().front();
        out.rows = allocate<std::optional<Row>>(memory, rows.size());
        auto moved = allocate<std::uint8_t>(memory, moving.size());
        a::LocalRootScratch scratch;
        for (std::size_t i = 0; i < rows.size(); ++i) {
            work.step();
            const auto& row = rows[i];
            const auto facet = row.material_point.facet;
            const auto panel = take(contact.uvlm_quad(facet, work.limit - work.used, stop));
            charge(work, panel.used_steps);
            if (panel.role != a::AeroFaceRole::Bound)
                reject(InvalidInput);
            for (const auto child : panel.facets) {
                for (const auto vertex : take(geometry.facet_points(child))) {
                    const auto column = moving_index(moving, take(geometry.owner(vertex)), work);
                    if (column == moving.size())
                        reject(InvalidInput);
                    moved.values()[column] = 1;
                }
            }
            auto material = resolve(geometry, {row.query, row.material_point});
            const auto xyz = position(geometry, material, a::Configuration::Current);
            const auto motion_index = moving.size() + i;
            const auto velocity = [&](std::size_t k) {
                return centre(take(motion.point(motion_index, k, p::JetOrder::First)));
            };
            for (std::size_t k = 0; k < axes; ++k) {
                const auto location =
                    m::exact(xyz.numerator[k],
                             -static_cast<int>(a::ExactQueryCoordinates::denominator_shift));
                qualify_point(location, take(motion.point(motion_index, k, p::JetOrder::Value)),
                              projection.point[k], limits.motion.maximum_interval_input_bits, work);
            }
            out.rows.values()[i].emplace(Row{std::move(material),
                                             tangents(geometry, facet, scratch, work),
                                             {velocity(0), velocity(1), velocity(2)},
                                             row});
        }
        work.step();
        for (const auto checked : moved.values()) {
            work.step();
            if (checked == 0)
                reject(InvalidInput);
        }
        std::ranges::sort(point_rows.values(), [&](const auto& left, const auto& right) {
            work.step();
            return left.id < right.id;
        });
        out.points = std::move(point_rows);
        out.moving_count = moving.size();
        work.step();
        out.contact.emplace(std::move(contact));
        out.motion.emplace(std::move(motion));
        return PrescribedAeroTrial{std::move(owner)};
    }
};
} // namespace aero_attachment_detail

PrescribedAeroTrial::PrescribedAeroTrial(
    core::PlannedBuffer<aero_attachment_detail::Trial>&& trial) noexcept
    : trial_(std::move(trial)) {}
PrescribedAeroTrial::PrescribedAeroTrial(PrescribedAeroTrial&&) noexcept = default;
PrescribedAeroTrial& PrescribedAeroTrial::operator=(PrescribedAeroTrial&&) noexcept = default;
PrescribedAeroTrial::~PrescribedAeroTrial() = default;

core::Result<std::reference_wrapper<const CaseMotionEvaluation>>
PrescribedAeroTrial::motion() const& noexcept {
    const auto data = trial_.values();
    if (!data.empty()) {
        if (const auto& motion = data.front().motion)
            return std::cref(*motion);
    }
    return std::unexpected(core::Failure{InvalidInput, evaluation});
}
core::Result<std::reference_wrapper<const a::GlobalContactState>>
PrescribedAeroTrial::contact() const& noexcept {
    const auto data = trial_.values();
    if (!data.empty()) {
        if (const auto& contact = data.front().contact)
            return std::cref(*contact);
    }
    return std::unexpected(core::Failure{InvalidInput, evaluation});
}

core::Result<io::SiEnclosure>
PrescribedAeroTrial::material_point(a::PointId id, std::size_t axis,
                                    p::JetOrder derivative) const noexcept {
    const auto data = trial_.values();
    if (!data.empty()) {
        const auto& run = data.front();
        const auto points = run.points.values();
        const auto row =
            std::ranges::lower_bound(points, id, {}, &aero_attachment_detail::PointOrdinal::id);
        if (row != points.end() && row->id == id && row->column < run.moving_count) {
            if (const auto& sampled = run.motion)
                return sampled->point(row->column, axis, derivative);
        }
    }
    return std::unexpected(core::Failure{InvalidInput, evaluation});
}

std::size_t PrescribedAeroTrial::collocation_count() const noexcept {
    const auto data = trial_.values();
    return data.empty() ? 0 : data.front().rows.size();
}
core::Result<DiscreteCollocation> PrescribedAeroTrial::collocation(std::size_t i) const& noexcept {
    const auto data = trial_.values();
    if (!data.empty()) {
        const auto rows = data.front().rows.values();
        if (i < rows.size()) {
            if (const auto& row = rows[i])
                return DiscreteCollocation{row->source.material_point.facet,
                                           row->point,
                                           row->surface,
                                           row->velocity,
                                           row->source.circulation,
                                           row->source.circulation_direction};
        }
    }
    return std::unexpected(core::Failure{InvalidInput, evaluation});
}

core::Result<PrescribedAeroAttachment> prepare_prescribed_aero_attachment(
    const io::SiCase& source, const io::schema::AeroAttachmentPrescribed& attachment,
    CaseMotionTimeWindow window, p::PchipInputStamp stamp, p::CurveWorkspace& curves,
    const CaseMotionLimits& limits, core::ResourcePlan& memory, const std::stop_token& stop,
    core::FailureDetails details) noexcept {
    return core::guard_entry(core::Operation::PrepareModel, PrecisionUnresolved, details,
                             [&]() -> core::Result<PrescribedAeroAttachment> {
                                 try {
                                     return aero_attachment_detail::Builder::prepare(
                                         source, attachment, window, stamp, curves, limits, memory,
                                         stop, details);
                                 } catch (const core::Failure& failure) {
                                     return std::unexpected(failure);
                                 }
                             });
}

core::Result<PrescribedAeroTrial> evaluate_prescribed_aero_attachment(
    const PrescribedAeroAttachment& prepared, const PrescribedAeroAttachment& current,
    const PrescribedAeroReference& reference, CaseMotionTimeQuery time,
    const PrescribedAeroProjection& projection, p::CurveWorkspace& curves,
    a::GeometryWorkspace& coordinates, a::GlobalWorkspace& global,
    const PrescribedAeroLimits& limits, core::ResourcePlan& memory, const std::stop_token& stop,
    core::FailureDetails details) noexcept {
    return core::guard_entry(
        evaluation, PrecisionUnresolved, details, [&]() -> core::Result<PrescribedAeroTrial> {
            try {
                return aero_attachment_detail::Builder::evaluate(
                    prepared, current, reference, time, projection, curves, coordinates, global,
                    limits, memory, stop, details);
            } catch (const core::Failure& failure) {
                return std::unexpected(failure);
            } catch (const m::ArithmeticFailure& failure) {
                return std::unexpected(core::Failure{failure.code, evaluation, details});
            }
        });
}

} // namespace uvlas::coupling::detail
