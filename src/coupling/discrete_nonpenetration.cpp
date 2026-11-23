#include "discrete_nonpenetration.hpp"

#include "aero/local_geometry_bounds_private.hpp"
#include "aero/p1_interval_private.hpp"
#include "aero/p1_shared_kernel.hpp"

#include <algorithm>
#include <bit>
#include <limits>

namespace uvlas::coupling::detail {
namespace a = aero::detail;
namespace p = a::p1;
namespace {
constexpr auto evaluate = core::Operation::EvaluateAero;
using Vector = std::array<p::Interval, 3>;

[[noreturn]] void reject(core::FailureCode code) {
    throw p::ArithmeticFailure{code};
}
template <class Result>
auto require(Result result) {
    if (!result) {
        const auto& error = result.error();
        if constexpr (std::is_same_v<typename Result::error_type, core::Failure>)
            throw error;
        else if constexpr (std::is_same_v<typename Result::error_type, a::GlobalGeometryFailure>)
            throw error.origin.cause;
        else
            throw error.cause;
    }
    return std::move(*result);
}
std::uint64_t remaining(const p::Work& work) noexcept {
    return work.limit - work.used;
}
std::uint64_t sum_budget(std::uint64_t a, std::uint64_t b) {
    if (b > std::numeric_limits<std::uint64_t>::max() - a)
        reject(core::FailureCode::ResourceExceeded);
    return a + b;
}
void account(p::Work& work, std::uint64_t amount) {
    if (amount > remaining(work))
        reject(core::FailureCode::ResourceExceeded);
    work.used += amount;
}
template <class T>
core::PlannedBuffer<T> allocate(core::ResourcePlan& resources, std::uint64_t count) {
    const auto size = static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max());
    if (count > size)
        reject(core::FailureCode::ResourceExceeded);
    return require(resources.allocate<T>(static_cast<core::Index>(count), evaluate));
}

p::Interval import(const a::ScalarWords& words, p::Work& work) {
    work.step();
    if (const auto* single = std::get_if<a::SingleWord>(&words))
        return a::decode_p1_word(single->value);
    const auto* two = std::get_if<a::DoubleWord>(&words);
    if (two == nullptr)
        reject(core::FailureCode::InvalidInput);
    return p::add(a::decode_p1_word(two->high), a::decode_p1_word(two->low), work);
}
Vector import(const NonpenetrationVector& value, p::Work& work) {
    return {import(value[0], work), import(value[1], work), import(value[2], work)};
}
p::Interval enclose(p::Interval center, double radius, int scale, p::Work& work) {
    const auto error = a::decode_p1_word(radius);
    return p::scale(
        p::Interval{p::subtract(center, error, work).lower, p::add(center, error, work).upper},
        scale);
}
p::Interval import(const a::ScaledEnclosure& value, p::Work& work) {
    return enclose(import(value.center(), work), value.radius(), value.exponent(), work);
}
p::Interval import(const io::SiEnclosure& value, p::Work& work) {
    return std::visit(
        [&](const auto& word) {
            using S = std::decay_t<decltype(word)>;
            const auto center = [&]() -> a::ScalarWords {
                if constexpr (std::is_same_v<S, io::ScaledBinary64>)
                    return a::SingleWord{word.center};
                else
                    return a::DoubleWord{word.high, word.low};
            }();
            return enclose(import(center, work), word.radius, p::exponent(word.exponent2), work);
        },
        value);
}
Vector cross(const Vector& a, const Vector& b, p::Work& work) {
    Vector result{};
    for (std::size_t axis = 0; axis < result.size(); ++axis) {
        const auto next = (axis + 1) % 3;
        const auto last = (axis + 2) % 3;
        result[axis] = p::subtract(p::multiply(a[next], b[last], work),
                                   p::multiply(a[last], b[next], work), work);
    }
    return result;
}
p::Interval dot(const Vector& a, const Vector& b, p::Work& work) {
    p::Interval result{};
    for (std::size_t axis = 0; axis < a.size(); ++axis)
        result = p::add(result, p::multiply(a[axis], b[axis], work), work);
    return result;
}
Vector difference(const Vector& a, const Vector& b, p::Work& work) {
    return {p::subtract(a[0], b[0], work), p::subtract(a[1], b[1], work),
            p::subtract(a[2], b[2], work)};
}
Vector qualify_normal(const a::OwnedGeometry& geometry, const a::UvlmQuad& panel,
                      const DiscreteCollocation& input, p::Work& work) {
    const auto tx = import(input.tangents[0], work);
    const auto ty = import(input.tangents[1], work);
    auto normal = cross(tx, ty, work);
    const auto length2 = dot(normal, normal, work);
    if (p::compare(length2.lower, {}) <= 0)
        reject(core::FailureCode::GeometryInvalid);
    for (auto facet : panel.facets) {
        const auto ids = require(geometry.facet_points(facet));
        const auto x0 =
            import(require(geometry.point_words(ids[0], a::Configuration::Current)), work);
        const auto x1 =
            import(require(geometry.point_words(ids[1], a::Configuration::Current)), work);
        const auto x2 =
            import(require(geometry.point_words(ids[2], a::Configuration::Current)), work);
        const auto edge1 = difference(x1, x0, work);
        const auto edge2 = difference(x2, x0, work);
        const auto outward = cross(edge1, edge2, work);
        if (p::compare(dot(normal, outward, work).lower, {}) <= 0)
            reject(core::FailureCode::GeometryInvalid);
    }
    const auto length = p::square_root(length2, work);
    for (auto& coordinate : normal)
        coordinate = p::divide(coordinate, length, work);
    return normal;
}

// The geometry context presently holds a binary64 time. Compare the exact
// dyadic, using the canonical limbs of the independently qualified flow jet.
bool same_time(double word, a::pchip::RationalView time, p::Work& work) {
    work.step();
    const auto value = a::decode_p1_word(word).lower;
    if (value.mantissa == 0)
        return time.sign == 0;
    const bool negative = value.mantissa < 0;
    if (time.sign != (negative ? -1 : 1))
        return false;
    const auto magnitude = p::magnitude(value.mantissa);
    const auto trailing = boost::multiprecision::lsb(magnitude);
    const p::Integer normalized = magnitude >> trailing;
    const auto mantissa = normalized.convert_to<std::uint64_t>();
    const auto power = p::exponent(static_cast<std::int64_t>(value.power) + trailing);
    const auto equal_words = [&](std::span<const std::uint32_t> limbs, std::uint64_t digits,
                                 unsigned shift) {
        const auto width = static_cast<unsigned>(std::bit_width(digits)) + shift;
        if (limbs.size() != (width + 31U) / 32U)
            return false;
        const auto first = shift / 32U;
        const auto rem = shift % 32U;
        const std::array<std::uint32_t, 3> words{
            static_cast<std::uint32_t>(digits << rem),
            static_cast<std::uint32_t>((digits << rem) >> 32U),
            rem == 0 ? 0U : static_cast<std::uint32_t>(digits >> (64U - rem))};
        for (std::size_t i = 0; i < limbs.size(); ++i) {
            work.step();
            if (limbs[i] != (i < first ? 0U : words[i - first]))
                return false;
        }
        return true;
    };
    return equal_words(time.numerator, mantissa, static_cast<unsigned>(std::max(power, 0))) &&
           equal_words(time.denominator, 1, static_cast<unsigned>(std::max(-power, 0)));
}
struct PanelRow {
    a::UvlmQuad panel{
        a::PhysicalFaceId{0}, a::AeroFaceRole::Bound, {a::FacetId{0}, a::FacetId{0}}, {}, 0};
    std::size_t caller = 0;
};
struct Trace {
    a::P1StrengthId id{0};
    a::ScalarWords value{a::SingleWord{0}};
    a::ScalarWords direction{a::SingleWord{0}};
};
a::ScalarWords negative(a::ScalarWords words) {
    if (const auto* one = std::get_if<a::SingleWord>(&words))
        return a::SingleWord{-one->value};
    const auto* two = std::get_if<a::DoubleWord>(&words);
    if (two == nullptr)
        reject(core::FailureCode::InvalidInput);
    return a::DoubleWord{-two->high, -two->low};
}
template <class T, class Key, class Projection>
const T& find(std::span<const T> items, Key key, Projection project, p::Work& work) {
    const auto it =
        std::lower_bound(items.begin(), items.end(), key, [&](const T& value, Key query) {
            work.step();
            return std::invoke(project, value) < query;
        });
    if (it == items.end() || std::invoke(project, *it) != key)
        reject(core::FailureCode::InvalidInput);
    return *it;
}
} // namespace

struct NonpenetrationStorage {
    std::optional<DiscreteNonpenetrationRow> value;
};
DiscreteNonpenetration::DiscreteNonpenetration(core::PlannedBuffer<NonpenetrationStorage> storage,
                                               ProjectedCaseUniformFlow&& flow,
                                               std::optional<a::P1ActionValue> geometry_owner,
                                               a::GeometryContext geometry_context,
                                               a::P1ActionStamp trial, std::uint64_t steps) noexcept
    : rows_(std::move(storage)),
      flow_(std::move(flow)),
      geometry_owner_(std::move(geometry_owner)),
      geometry_context_(geometry_context),
      trial_stamp_(trial),
      used_steps_(steps) {}
DiscreteNonpenetration::DiscreteNonpenetration(DiscreteNonpenetration&&) noexcept = default;
DiscreteNonpenetration&
DiscreteNonpenetration::operator=(DiscreteNonpenetration&&) noexcept = default;
DiscreteNonpenetration::~DiscreteNonpenetration() = default;
std::size_t DiscreteNonpenetration::size() const noexcept {
    return rows_.size();
}
core::Result<DiscreteNonpenetrationRow> DiscreteNonpenetration::row(std::size_t i) const noexcept {
    if (i >= rows_.size())
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, evaluate});
    const auto& result = rows_.values()[i].value;
    if (!result)
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, evaluate});
    return *result;
}
const ProjectedCaseUniformFlow& DiscreteNonpenetration::flow() const& noexcept {
    return flow_;
}
std::uint64_t DiscreteNonpenetration::used_steps() const noexcept {
    return flow_.source().jet.active() && (rows_.size() == 0 || geometry_owner_) ? used_steps_ : 0;
}

a::GeometryContext DiscreteNonpenetration::geometry_context() const noexcept {
    return geometry_context_;
}
a::P1ActionStamp DiscreteNonpenetration::trial_stamp() const noexcept {
    return trial_stamp_;
}

struct NonpenetrationBuilder {
    static DiscreteNonpenetration
    evaluate_rows(const DiscreteNonpenetrationInput& input, ProjectedCaseUniformFlow&& flow,
                  const DiscreteNonpenetrationLimits& limits, a::SharedP1Workspace& shared,
                  a::LocalPacketWorkspace& local, a::P1Workspace& kernel,
                  core::ResourcePlan& memory, const std::stop_token& stop) {
        p::Work work{0, limits.steps, stop};
        work.step();
        if (!a::p1_environment_supported())
            reject(core::FailureCode::PrecisionUnsupported);
        if (const auto invalid = limits.geometry_accuracy.validation_error())
            reject(*invalid);
        if (const auto invalid = limits.induction_accuracy.validation_error())
            reject(*invalid);
        const auto geometry_owner = require(input.contact.geometry());
        const auto& geometry = geometry_owner.get();
        const auto& jet = flow.source().jet;
        if (require(jet.stamp()) != input.flow_stamp ||
            flow.source().requested_side != input.flow_side ||
            !same_time(geometry.context().time, require(jet.time()), work))
            reject(core::FailureCode::StaleGeneration);
        if (!input.wake.empty() && geometry.context().diagnostics.parent != input.history.parent)
            reject(core::FailureCode::StaleGeneration);

        const auto counts = require(geometry.counts());
        const auto face_count = counts.faces;
        if (counts.quads != face_count || counts.triangles != 0)
            reject(core::FailureCode::InvalidInput);
        const auto bound_count = input.collocations.size();
        if (bound_count > face_count || input.wake.size() != face_count - bound_count ||
            input.contact.trailing_edge_count() > input.wake.size())
            reject(core::FailureCode::InvalidInput);
        const auto maximum = static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max());
        if (face_count > maximum || input.history.columns.size() > maximum - face_count)
            reject(core::FailureCode::ResourceExceeded);
        auto rows = allocate<NonpenetrationStorage>(memory, bound_count);
        auto bound = allocate<PanelRow>(memory, bound_count);
        auto maps = allocate<a::P1FaceStrengths>(memory, face_count);
        auto traces = allocate<Trace>(memory, face_count + input.history.columns.size());
        auto attachments = allocate<std::uint64_t>(memory, input.contact.trailing_edge_count());
        auto trace = traces.values();
        std::size_t num_traces = 0;
        for (std::size_t i = 0; i < bound_count; ++i) {
            const auto& row = input.collocations[i];
            const auto panel = require(input.contact.uvlm_quad(row.facet, remaining(work), stop));
            account(work, panel.used_steps);
            if (panel.role != a::AeroFaceRole::Bound)
                reject(core::FailureCode::InvalidInput);
            bound.values()[i] = PanelRow{panel, i};
            maps.values()[i] = a::P1FaceStrengths{panel.face, a::P0BoundColumn{}};
            trace[num_traces++] = {a::P1StrengthId{panel.face.value()}, row.circulation,
                                   row.circulation_direction};
        }
        std::ranges::sort(bound.values(), [&](const PanelRow& l, const PanelRow& r) {
            work.step();
            return l.panel.face < r.panel.face;
        });
        auto prior = input.history.columns.end();
        for (auto it = input.history.columns.begin(); it != input.history.columns.end(); ++it) {
            work.step();
            if (prior != input.history.columns.end() && prior->id >= it->id)
                reject(core::FailureCode::InvalidInput);
            prior = it;
            trace[num_traces++] = {it->id, it->strength, a::SingleWord{0}};
        }
        const auto order = [&](const Trace& l, const Trace& r) {
            work.step();
            return l.id < r.id;
        };
        std::ranges::sort(trace.first(num_traces), order);
        for (std::size_t i = 1; i < num_traces; ++i) {
            work.step();
            if (trace[i - 1].id == trace[i].id)
                reject(core::FailureCode::InvalidInput);
        }
        const auto existing = num_traces;
        std::size_t visit = 0, attached = 0;
        std::uint64_t spare = 0;
        const auto new_column = [&] {
            for (; visit < existing && trace[visit].id.value() <= spare; ++visit) {
                work.step();
                if (trace[visit].id.value() == spare)
                    ++spare;
            }
            return a::P1StrengthId{spare++};
        };
        const auto history = [&](a::P1StrengthId id) {
            return find(input.history.columns, id, &FrozenWakeTrace::id, work).id;
        };
        for (std::size_t i = 0; i < input.wake.size(); ++i) {
            const auto& row = input.wake[i];
            const auto panel = require(input.contact.uvlm_quad(row.facet, remaining(work), stop));
            account(work, panel.used_steps);
            const auto* previous = std::get_if<a::P1StrengthId>(&row.leading);
            if (panel.role != a::AeroFaceRole::Wake ||
                panel.attachment.has_value() == (previous != nullptr))
                reject(core::FailureCode::InvalidInput);
            auto leading = a::P1StrengthId{0};
            if (previous != nullptr) {
                leading = history(*previous);
            } else {
                if (!panel.attachment || attached >= attachments.size())
                    reject(core::FailureCode::InvalidInput);
                const auto& entry = *panel.attachment;
                attachments.values()[attached++] = entry.edge.value();
                const auto& source = find(
                    std::span<const PanelRow>{bound.values()}, entry.bound_face,
                    [](const PanelRow& b) { return b.panel.face; }, work);
                leading = a::P1StrengthId{entry.bound_face.value()};
                if (entry.coefficient == -1) {
                    const auto& collocation = input.collocations[source.caller];
                    leading = new_column();
                    trace[num_traces++] = {leading, negative(collocation.circulation),
                                           negative(collocation.circulation_direction)};
                }
            }
            maps.values()[bound_count + i] =
                a::P1FaceStrengths{panel.face, std::array{leading, history(row.trailing)}};
        }
        if (attached != attachments.size())
            reject(core::FailureCode::InvalidInput);
        std::ranges::sort(attachments.values(), [&](std::uint64_t l, std::uint64_t r) {
            work.step();
            return l < r;
        });
        for (std::size_t i = 1; i < attachments.size(); ++i) {
            work.step();
            if (attachments.values()[i - 1] == attachments.values()[i])
                reject(core::FailureCode::InvalidInput);
        }
        std::ranges::sort(maps.values(),
                          [&](const a::P1FaceStrengths& l, const a::P1FaceStrengths& r) {
                              work.step();
                              return l.face < r.face;
                          });
        for (std::size_t i = 1; i < maps.size(); ++i) {
            work.step();
            if (maps.values()[i - 1].face == maps.values()[i].face)
                reject(core::FailureCode::InvalidInput);
        }
        std::ranges::sort(trace.first(num_traces), order);
        auto ids = allocate<a::P1StrengthId>(memory, num_traces);
        auto values = allocate<a::P1ComplexWords>(memory, num_traces);
        auto directions = allocate<a::P1ComplexWords>(memory, num_traces);
        for (std::size_t i = 0; i < num_traces; ++i) {
            work.step();
            (void)import(trace[i].value, work);
            (void)import(trace[i].direction, work);
            ids.values()[i] = trace[i].id;
            values.values()[i].real = trace[i].value;
            directions.values()[i].real = trace[i].direction;
        }
        std::uint64_t payload = 0, reserved_steps = 0;
        const auto application_reserve = sum_budget(limits.action_steps, limits.action_steps);
        for (const auto& item : bound.values()) {
            const auto& row = input.collocations[item.caller];
            const a::P1ActionSource source{input.contact, row.query, input.trial_stamp};
            auto preflight_limits = limits.induction;
            preflight_limits.steps = std::min(preflight_limits.steps, remaining(work));
            const auto plan = require(a::forecast_bound_wake_induction(
                source, maps.values(), ids.values(), preflight_limits, stop));
            account(work, plan.preflight_steps);
            payload = sum_budget(payload, plan.payload_bytes);
            reserved_steps =
                sum_budget(reserved_steps, sum_budget(plan.reserved_steps, application_reserve));
        }
        if (payload > memory.remaining_bytes() || reserved_steps > remaining(work))
            reject(core::FailureCode::ResourceExceeded);

        const Vector ambient{import(require(flow.velocity(0)), work),
                             import(require(flow.velocity(1)), work),
                             import(require(flow.velocity(2)), work)};
        std::optional<a::P1ActionValue> geometry_retainer;
        for (std::size_t i = 0; i < bound_count; ++i) {
            const auto& item = bound.values()[i];
            const auto& row = input.collocations[item.caller];
            const auto normal = qualify_normal(geometry, item.panel, row, work);
            const auto relative = difference(ambient, import(row.velocity, work), work);
            const a::P1ActionSource source{input.contact, row.query, input.trial_stamp};
            auto induction_limits = limits.induction;
            induction_limits.steps = std::min(induction_limits.steps, remaining(work));
            const auto action = require(a::prepare_bound_wake_induction(
                source, maps.values(), ids.values(), limits.geometry_accuracy,
                limits.induction_accuracy, induction_limits, shared, local, kernel, memory, stop));
            const auto report = require(action.report());
            account(work, report.used_steps);
            auto full = require(a::apply_shared_p1(action, source, values.values(), std::nullopt,
                                                   std::min(limits.action_steps, remaining(work)),
                                                   shared, stop));
            account(work, full.used_steps());
            const auto linear = require(
                a::apply_shared_p1(action, source, directions.values(), std::nullopt,
                                   std::min(limits.action_steps, remaining(work)), shared, stop));
            account(work, linear.used_steps());
            Vector velocity{}, delta{};
            for (std::size_t axis = 0; axis < normal.size(); ++axis) {
                velocity[axis] =
                    p::add(import(full.velocity()[axis].real, work), relative[axis], work);
                delta[axis] = import(linear.velocity()[axis].real, work);
            }
            a::LocalRootScratch rounding;
            work.step();
            const auto value = a::enclose_p1_value(
                dot(normal, velocity, work), limits.induction_accuracy.representation, rounding);
            work.step();
            const auto direction = a::enclose_p1_value(
                dot(normal, delta, work), limits.induction_accuracy.representation, rounding);
            rows.values()[i].value.emplace(DiscreteNonpenetrationRow{
                item.panel.face, value, direction, report.far_facets != 0});
            if (i == 0)
                geometry_retainer.emplace(std::move(full));
        }
        work.step();
        return {std::move(rows),    std::move(flow),   std::move(geometry_retainer),
                geometry.context(), input.trial_stamp, work.used};
    }
};

core::Result<DiscreteNonpenetration> evaluate_discrete_nonpenetration(
    const DiscreteNonpenetrationInput& input, ProjectedCaseUniformFlow&& flow,
    const DiscreteNonpenetrationLimits& limits, a::SharedP1Workspace& shared,
    a::LocalPacketWorkspace& local, a::P1Workspace& kernel, core::ResourcePlan& memory,
    const std::stop_token& stop) noexcept {
    try {
        return NonpenetrationBuilder::evaluate_rows(input, std::move(flow), limits, shared, local,
                                                    kernel, memory, stop);
    } catch (const p::ArithmeticFailure& error) {
        return std::unexpected(core::Failure{error.code, evaluate});
    } catch (const core::Failure& error) {
        return std::unexpected(error);
    } catch (const std::bad_alloc&) {
        return std::unexpected(core::Failure{core::FailureCode::ResourceExceeded, evaluate});
    } catch (...) {
        return std::unexpected(core::Failure{core::FailureCode::PrecisionUnresolved, evaluate});
    }
}
} // namespace uvlas::coupling::detail
