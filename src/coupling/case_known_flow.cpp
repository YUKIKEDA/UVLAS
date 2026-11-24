#include "case_known_flow.hpp"

#include "aero/certified_trigonometry_private.hpp"
#include "aero/local_geometry_bounds_private.hpp"
#include "aero/p1_shared_kernel.hpp"
#include "case_component_precision_private.hpp"
#include "case_flow_series_private.hpp"
#include "io/si_exact_words.hpp"

#include <algorithm>

namespace uvlas::coupling::detail {
namespace a = aero::detail;
namespace p = a::pchip;
namespace m = a::p1;
namespace s = io::schema;
namespace {
using enum core::FailureCode;
constexpr auto prepare = core::Operation::PrepareModel;
constexpr auto evaluate = core::Operation::EvaluateAero;
constexpr std::size_t xyz = 3;
constexpr io::SchemaDimension unitless{{0, 0, 0, 0}}, metres{{1, 0, 0, 0}}, seconds{{0, 0, 1, 0}},
    wavevector{{-1, 0, 0, 1}}, frequency{{0, 0, -1, 1}}, angle{{0, 0, 0, 1}};

template <class T>
T take(core::Result<T> result) {
    if (!result)
        throw result.error();
    return std::move(*result);
}

template <class T, class F>
core::Result<T> guarded(core::Operation op, core::FailureDetails details, F call) noexcept {
    return core::guard_entry(op, PrecisionUnresolved, details, [&]() -> core::Result<T> {
        try {
            return call();
        } catch (const core::Failure& failure) {
            return std::unexpected(failure);
        } catch (const m::ArithmeticFailure& failure) {
            return std::unexpected(core::Failure{failure.code, op, details});
        }
    });
}

[[noreturn]] void reject(core::FailureCode code) {
    throw m::ArithmeticFailure{code};
}

struct Number {
    core::PlannedBuffer<std::uint32_t> words;
    std::size_t numerator = 0;
    int sign = 0;
    [[nodiscard]] p::RationalView view() const noexcept {
        return {sign, words.values().first(numerator), words.values().subspan(numerator)};
    }
};

Number own(p::RationalView value, core::ResourcePlan& memory, m::Work& work,
           core::Operation operation = prepare) {
    work.step();
    const auto size = value.numerator.size() + value.denominator.size();
    auto words = take(memory.allocate<std::uint32_t>(static_cast<core::Index>(size), operation));
    std::size_t i = 0;
    for (const auto part : {value.numerator, value.denominator})
        for (const auto limb : part) {
            work.step();
            words.values()[i++] = limb;
        }
    return {std::move(words), value.numerator.size(), value.sign};
}

using case_components::rational;
using case_components::validate_words;

constexpr std::size_t direction_at = 0, wave_at = 3, origin_at = 6, omega_at = 9, phase_at = 10,
                      time_origin_at = 11;
struct Gust {
    std::array<Number, 12> input;
    m::PhaseUnits units = m::PhaseUnits::Radians;
};

struct PhysicalGust {
    std::array<m::Interval, xyz> direction, wave, origin, wave_radians;
    m::Interval omega, phase, time_origin, omega_radians, pi;
};

std::optional<PhysicalGust> qualify(const std::optional<Gust>& gust, std::uint64_t input_bits,
                                    m::Work& work) {
    if (!gust)
        return std::nullopt;
    const auto at = [&](std::size_t i) {
        return rational(gust->input[i].view(), input_bits, work);
    };
    PhysicalGust result{};
    result.pi = m::pi_constant(work);
    const auto radians = [&](const m::Interval& argument) {
        return gust->units == m::PhaseUnits::PiMultiple ? m::multiply(argument, result.pi, work)
                                                        : argument;
    };
    result.omega = at(omega_at);
    result.omega_radians = radians(result.omega);
    result.phase = at(phase_at);
    result.time_origin = at(time_origin_at);
    auto norm_squared = m::integer(0);
    for (std::size_t i = 0; i < xyz; ++i) {
        result.direction[i] = at(direction_at + i);
        result.wave[i] = at(wave_at + i);
        result.wave_radians[i] = radians(result.wave[i]);
        result.origin[i] = at(origin_at + i);
        const auto& component = result.direction[i];
        norm_squared = m::add(norm_squared, m::multiply(component, component, work), work);
    }
    const auto norm = m::square_root(norm_squared, work);
    for (auto& component : result.direction)
        component = m::divide(component, norm, work);
    return result;
}

std::optional<Gust> clone(const std::optional<Gust>& source, std::uint64_t input_bits,
                          core::ResourcePlan& memory, m::Work& work) {
    if (!source)
        return std::nullopt;
    Gust saved{};
    saved.units = source->units;
    for (std::size_t i = 0; i < saved.input.size(); ++i) {
        const auto value = source->input[i].view();
        validate_words(value, input_bits);
        saved.input[i] = own(value, memory, work);
    }
    return saved;
}

void match(const std::optional<Gust>& field, const std::optional<Gust>& current,
           std::uint64_t input_bits, m::Work& work) {
    work.step();
    if (field.has_value() != current.has_value())
        reject(StaleGeneration);
    if (!field || !current)
        return;
    if (field->units != current->units)
        reject(StaleGeneration);
    const auto equal = [&](auto a_words, auto b_words) {
        return std::equal(a_words.begin(), a_words.end(), b_words.begin(), b_words.end(),
                          [&](auto l, auto r) {
                              work.step();
                              return l == r;
                          });
    };
    for (std::size_t i = 0; i < field->input.size(); ++i) {
        work.step();
        const auto l = field->input[i].view();
        const auto r = current->input[i].view();
        validate_words(l, input_bits);
        validate_words(r, input_bits);
        if (l.sign != r.sign || !equal(l.numerator, r.numerator) ||
            !equal(l.denominator, r.denominator))
            reject(StaleGeneration);
    }
}

std::optional<Gust> prepare_gust(const io::SiCase& source, const CaseKnownFlowLimits& limits,
                                 core::ResourcePlan& memory, m::Work& work) {
    const auto schema = take(source.source_case());
    const auto root = take(schema.get().root());
    const auto flow = root.flow.record();
    const auto* known = std::get_if<s::FlowKnown>(&flow);
    if (known == nullptr)
        reject(Unsupported);
    const auto selected = known->gust.record();
    const auto* sine = std::get_if<s::GustSine>(&selected);
    if (sine == nullptr)
        return std::nullopt;
    const auto& exact = limits.time_curves.conversion;
    const auto phase_recipe = take(source.quantity(sine->phase)).recipe;
    const bool degrees = phase_recipe.units.angle == io::SiAngleUnit::Degree;
    const auto d180 = degrees ? take(structure::ExactRational::from_decimal("180", exact, memory))
                              : structure::ExactRational{};
    const auto read = [&](auto field, io::SchemaDimension expected) {
        work.step();
        auto recipe = take(source.quantity(field)).recipe;
        if (recipe.dimensions.exponents != expected.exponents)
            reject(InvalidInput);
        const bool scale_angle = degrees && expected.exponents[3] == 1;
        if (scale_angle) {
            // Remove the common pi multiplier, retaining the exact coefficient of pi.
            recipe.units.angle = io::SiAngleUnit::Radian;
            recipe.pi_power = 0;
        }
        auto value = take(io::materialize_rational(recipe, exact, memory));
        if (scale_angle)
            return take(value.divide(d180, exact, memory));
        return value;
    };
    std::array<structure::ExactRational, 12> numbers{};
    structure::ExactRational orthogonality;
    bool nonzero_direction = false;
    for (std::size_t i = 0; i < xyz; ++i) {
        const auto at = static_cast<core::Index>(i);
        auto e = read(take(sine->direction.at(at)), unitless);
        auto k = read(take(sine->wavevector.at(at)), wavevector);
        nonzero_direction = nonzero_direction || e.sign() != 0;
        const auto term = take(e.multiply(k, exact, memory));
        orthogonality = take(orthogonality.add(term, exact, memory));
        numbers[direction_at + i] = std::move(e);
        numbers[wave_at + i] = std::move(k);
        numbers[origin_at + i] = read(take(sine->origin.at(at)), metres);
    }
    if (!nonzero_direction || orthogonality.sign() != 0)
        reject(InvalidInput);
    numbers[omega_at] = read(sine->angular_frequency, frequency);
    numbers[phase_at] = read(sine->phase, angle);
    numbers[time_origin_at] = read(sine->time_origin, seconds);
    Gust result{};
    result.units = degrees ? m::PhaseUnits::PiMultiple : m::PhaseUnits::Radians;
    for (std::size_t i = 0; i < numbers.size(); ++i) {
        if (numbers[i].numerator_bits() > limits.maximum_input_bits ||
            numbers[i].denominator_bits() > limits.maximum_input_bits ||
            numbers[i].numerator_bits() > 32768U || numbers[i].denominator_bits() > 32768U)
            reject(ResourceExceeded);
        const auto source_words =
            take(io::detail::copy_rational_words(numbers[i], memory, work.stop));
        const auto value = source_words.view();
        result.input[i] = own({value.sign, value.numerator, value.denominator}, memory, work);
    }
    return result;
}

using case_components::project;

void add_gust(const PhysicalGust& gust, m::PhaseUnits units,
              const a::ExactUniformFlowJet& amplitude, const CaseKnownFlowQuery& query,
              std::array<m::Interval, 15>& values, std::uint64_t input_bits, m::Work& work) {
    const auto theta_t =
        m::subtract(rational(query.time.time, input_bits, work), gust.time_origin, work);
    auto phase = m::subtract(gust.phase, m::multiply(gust.omega, theta_t, work), work);
    for (std::size_t i = 0; i < xyz; ++i) {
        const auto dx =
            m::subtract(rational(query.position[i], input_bits, work), gust.origin[i], work);
        phase = m::add(phase, m::multiply(gust.wave[i], dx, work), work);
    }
    const auto trig = m::sine_cosine(phase, units, gust.pi, work);
    const auto amount = rational(take(amplitude.velocity(0)), input_bits, work);
    const auto slope = rational(take(amplitude.time_partial(0)), input_bits, work);
    const auto gust_speed = m::multiply(amount, trig.sine, work);
    const auto ac = m::multiply(amount, trig.cosine, work);
    const auto gust_dt = m::subtract(m::multiply(slope, trig.sine, work),
                                     m::multiply(gust.omega_radians, ac, work), work);
    for (std::size_t i = 0; i < xyz; ++i) {
        const auto& e = gust.direction[i];
        values[i] = m::add(values[i], m::multiply(e, gust_speed, work), work);
        values[xyz + i] = m::add(values[xyz + i], m::multiply(e, gust_dt, work), work);
        const auto factor = m::multiply(e, ac, work);
        for (std::size_t j = 0; j < xyz; ++j)
            values[2 * xyz + xyz * i + j] = m::multiply(factor, gust.wave_radians[j], work);
    }
}
} // namespace
namespace known_flow_detail {

struct Prepared {
    std::optional<CaseUniformFlowReady> ambient, amplitude;
    std::optional<Gust> gust;
    std::optional<PhysicalGust> physical;
    std::optional<io::SchemaIdentity> free_wake_case;
};

struct Saved {
    std::optional<a::UniformFlowSnapshot> ambient, amplitude;
    std::optional<Gust> gust;
};

struct Position {
    std::array<Number, xyz> component;
};

struct Builder {
    static const Prepared& active(const CaseKnownFlowReady& field) {
        const auto records = field.data_.values();
        if (records.empty() || !records.front().ambient)
            reject(InvalidInput);
        return records.front();
    }

    static CaseKnownFlowReady
    publish(CaseUniformFlowReady ambient, std::optional<CaseUniformFlowReady> amplitude,
            std::optional<Gust> gust, std::optional<io::SchemaIdentity> free_wake_case,
            std::uint64_t input_bits, core::ResourcePlan& memory, m::Work& work) {
        auto physical = qualify(gust, input_bits, work);
        auto record = take(memory.allocate<Prepared>(1, prepare));
        record.values().front() = {std::move(ambient), std::move(amplitude), std::move(gust),
                                   std::move(physical), free_wake_case};
        work.step();
        return CaseKnownFlowReady{std::move(record)};
    }

    static CaseKnownFlowReady make(const io::SiCase& source, p::PchipInputStamp stamp,
                                   std::optional<CaseUniformTimeWindow> window,
                                   a::UniformFlowWorkspace& workspace,
                                   const CaseKnownFlowLimits& limits, core::ResourcePlan& memory,
                                   const std::stop_token& stop, core::FailureDetails details) {
        m::Work work{0, limits.interval_steps, stop};
        work.step();
        auto gust = prepare_gust(source, limits, memory, work);
        auto series = take(case_flow_detail::prepare_known_series(
            source, stamp, window, workspace, limits.time_curves, memory, stop, details));
        if (series.amplitude.has_value() != gust.has_value())
            reject(InvalidInput);
        const auto& schema = take(source.source_case()).get();
        const auto root = take(schema.root());
        const auto analysis = root.analysis.record();
        std::optional<io::SchemaIdentity> free_wake_case;
        if (const auto* transient = std::get_if<s::AnalysisTransient>(&analysis)) {
            const auto wake = transient->wake.record();
            if (const auto* active = std::get_if<s::WakeActive>(&wake);
                active != nullptr &&
                std::holds_alternative<s::ConvectionFree>(active->convection.record()))
                free_wake_case = take(schema.identity());
        }
        return publish(std::move(series.ambient), std::move(series.amplitude), std::move(gust),
                       free_wake_case, limits.maximum_input_bits, memory, work);
    }

    static CaseKnownFlowJet run(const CaseKnownFlowReady& requested, const CaseKnownFlowReady& host,
                                const CaseKnownFlowQuery& query,
                                const CaseKnownFlowProjection& projection,
                                a::UniformFlowWorkspace& workspace,
                                const CaseKnownFlowLimits& limits, core::ResourcePlan& memory,
                                const std::stop_token& stop, core::FailureDetails details) {
        m::Work work{0, limits.interval_steps, stop};
        work.step();
        if (!a::p1_environment_supported() ||
            (projection.representation != core::Precision::Binary64 &&
             projection.representation != core::Precision::TwoWord))
            reject(PrecisionUnsupported);

        const auto& field = active(requested);
        const auto& current = active(host);
        if (!field.ambient || !current.ambient)
            reject(InvalidInput);
        match(field.gust, current.gust, limits.maximum_input_bits, work);
        for (const auto& position : query.position) {
            validate_words(position, limits.maximum_input_bits);
            work.step();
        }
        const auto sample = [&](const CaseUniformFlowReady& saved,
                                const CaseUniformFlowReady& live) {
            return take(evaluate_case_uniform_flow(saved, live, query.time, workspace,
                                                   limits.time_curves.flow, memory, stop, details));
        };
        auto ambient = sample(*field.ambient, *current.ambient);
        std::array<m::Interval, 15> values{};
        for (std::size_t i = 0; i < xyz; ++i) {
            values[i] = rational(take(ambient.jet.velocity(i)), limits.maximum_input_bits, work);
            values[xyz + i] =
                rational(take(ambient.jet.time_partial(i)), limits.maximum_input_bits, work);
        }
        if (field.gust) {
            if (!field.amplitude || !current.amplitude || !field.physical)
                reject(InvalidInput);
            const auto amplitude = sample(*field.amplitude, *current.amplitude);
            add_gust(*field.physical, field.gust->units, amplitude.jet, query, values,
                     limits.maximum_input_bits, work);
        }
        a::LocalRootScratch scratch{};
        std::array<io::SiEnclosure, 15> projected{};
        const auto seal = [&](std::size_t i, const CaseFlowComponentBudget& budget) {
            projected[i] = project(values[i], budget, projection.representation,
                                   limits.maximum_input_bits, scratch, work);
        };
        for (std::size_t i = 0; i < xyz; ++i) {
            seal(i, projection.velocity[i]);
            seal(xyz + i, projection.time_partial[i]);
        }
        for (std::size_t i = 0; i < xyz * xyz; ++i)
            seal(2 * xyz + i, projection.spatial_partial[i]);
        auto position = take(memory.allocate<Position>(1, evaluate));
        for (std::size_t i = 0; i < xyz; ++i)
            position.values().front().component[i] = own(query.position[i], memory, work, evaluate);
        work.step();
        return CaseKnownFlowJet{std::move(ambient), projected, std::move(position), work.used};
    }

    static CaseKnownFlowSnapshot save(const CaseKnownFlowReady& requested,
                                      a::UniformFlowWorkspace& workspace,
                                      const CaseKnownFlowLimits& limits, core::ResourcePlan& memory,
                                      const std::stop_token& stop, core::FailureDetails details) {
        m::Work work{0, limits.interval_steps, stop};
        work.step();
        const auto& field = active(requested);
        if (!field.ambient)
            reject(InvalidInput);
        const auto capture = [&](const CaseUniformFlowReady& owner) {
            return take(snapshot_case_uniform_flow(owner, workspace, limits.time_curves.flow,
                                                   memory, stop, details));
        };
        auto ambient = capture(*field.ambient);
        std::optional<a::UniformFlowSnapshot> amplitude;
        if (field.amplitude)
            amplitude = capture(*field.amplitude);
        auto gust = clone(field.gust, limits.maximum_input_bits, memory, work);
        auto record = take(memory.allocate<Saved>(1, prepare));
        record.values().front() = {std::move(ambient), std::move(amplitude), std::move(gust)};
        work.step();
        return CaseKnownFlowSnapshot{std::move(record)};
    }

    static CaseKnownFlowReady restore(const CaseKnownFlowSnapshot& snapshot,
                                      const CaseKnownFlowReady& host,
                                      a::UniformFlowWorkspace& workspace,
                                      const CaseKnownFlowLimits& limits, core::ResourcePlan& memory,
                                      const std::stop_token& stop, core::FailureDetails details) {
        m::Work work{0, limits.interval_steps, stop};
        work.step();
        const auto& current = active(host);
        const auto records = snapshot.data_.values();
        if (records.empty())
            reject(InvalidInput);
        const auto& field = records.front();
        if (!field.ambient || !current.ambient)
            reject(InvalidInput);
        match(field.gust, current.gust, limits.maximum_input_bits, work);
        const auto requalify = [&](const a::UniformFlowSnapshot& saved,
                                   const CaseUniformFlowReady& live) {
            return take(restore_case_uniform_flow(take(saved.view()), live, workspace,
                                                  limits.time_curves.flow, memory, stop, details));
        };
        auto ambient = requalify(*field.ambient, *current.ambient);
        std::optional<CaseUniformFlowReady> amplitude;
        if (field.gust) {
            if (!field.amplitude || !current.amplitude)
                reject(InvalidInput);
            amplitude = requalify(*field.amplitude, *current.amplitude);
        }
        auto gust = clone(field.gust, limits.maximum_input_bits, memory, work);
        return publish(std::move(ambient), std::move(amplitude), std::move(gust),
                       current.free_wake_case, limits.maximum_input_bits, memory, work);
    }
};
} // namespace known_flow_detail

CaseKnownFlowReady::CaseKnownFlowReady(
    core::PlannedBuffer<known_flow_detail::Prepared>&& data) noexcept
    : data_(std::move(data)) {}
CaseKnownFlowReady::CaseKnownFlowReady(CaseKnownFlowReady&&) noexcept = default;
CaseKnownFlowReady& CaseKnownFlowReady::operator=(CaseKnownFlowReady&&) noexcept = default;
CaseKnownFlowReady::~CaseKnownFlowReady() = default;

core::Result<io::SchemaIdentity> CaseKnownFlowReady::free_wake_case_identity() const noexcept {
    const auto records = data_.values();
    if (records.empty())
        return std::unexpected(core::Failure{InvalidInput, evaluate});
    if (const auto& identity = records.front().free_wake_case)
        return *identity;
    return std::unexpected(core::Failure{Unsupported, evaluate});
}

CaseKnownFlowSnapshot::CaseKnownFlowSnapshot(
    core::PlannedBuffer<known_flow_detail::Saved>&& data) noexcept
    : data_(std::move(data)) {}
CaseKnownFlowSnapshot::CaseKnownFlowSnapshot(CaseKnownFlowSnapshot&&) noexcept = default;
CaseKnownFlowSnapshot& CaseKnownFlowSnapshot::operator=(CaseKnownFlowSnapshot&&) noexcept = default;
CaseKnownFlowSnapshot::~CaseKnownFlowSnapshot() = default;

CaseKnownFlowJet::CaseKnownFlowJet(CaseUniformFlowEvaluation&& time,
                                   const std::array<io::SiEnclosure, 15>& values,
                                   core::PlannedBuffer<known_flow_detail::Position>&& position,
                                   std::uint64_t steps) noexcept
    : time_(std::move(time)), position_(std::move(position)), values_(values), steps_(steps) {}
CaseKnownFlowJet::CaseKnownFlowJet(CaseKnownFlowJet&&) noexcept = default;
CaseKnownFlowJet& CaseKnownFlowJet::operator=(CaseKnownFlowJet&&) noexcept = default;
CaseKnownFlowJet::~CaseKnownFlowJet() = default;
core::Result<p::PchipInputStamp> CaseKnownFlowJet::stamp() const noexcept {
    return time_.jet.stamp();
}
core::Result<p::RationalView> CaseKnownFlowJet::position(std::size_t i) const& noexcept {
    if (position_.values().empty() || i >= xyz)
        return std::unexpected(core::Failure{InvalidInput, evaluate});
    return position_.values().front().component[i].view();
}
core::Result<p::RationalView> CaseKnownFlowJet::time() const& noexcept {
    return time_.jet.time();
}
p::RequestedSide CaseKnownFlowJet::requested_side() const noexcept {
    return time_.requested_side;
}
core::Result<io::SiEnclosure> CaseKnownFlowJet::velocity(std::size_t i) const noexcept {
    if (!time_.jet.active() || i >= xyz)
        return std::unexpected(core::Failure{InvalidInput, evaluate});
    return values_[i];
}
core::Result<io::SiEnclosure> CaseKnownFlowJet::time_partial(std::size_t i) const noexcept {
    if (!time_.jet.active() || i >= xyz)
        return std::unexpected(core::Failure{InvalidInput, evaluate});
    return values_[xyz + i];
}
core::Result<io::SiEnclosure> CaseKnownFlowJet::spatial_partial(std::size_t i,
                                                                std::size_t j) const noexcept {
    if (!time_.jet.active() || i >= xyz || j >= xyz)
        return std::unexpected(core::Failure{InvalidInput, evaluate});
    return values_[2 * xyz + xyz * i + j];
}

core::Result<CaseKnownFlowReady>
prepare_case_known_flow(const io::SiCase& source, p::PchipInputStamp stamp,
                        std::optional<CaseUniformTimeWindow> window,
                        a::UniformFlowWorkspace& workspace, const CaseKnownFlowLimits& limits,
                        core::ResourcePlan& memory, const std::stop_token& stop,
                        core::FailureDetails details) noexcept {
    return guarded<CaseKnownFlowReady>(prepare, details, [&] {
        return known_flow_detail::Builder::make(source, stamp, window, workspace, limits, memory,
                                                stop, details);
    });
}
core::Result<CaseKnownFlowJet>
evaluate_case_known_flow(const CaseKnownFlowReady& requested, const CaseKnownFlowReady& current,
                         const CaseKnownFlowQuery& query, const CaseKnownFlowProjection& projection,
                         a::UniformFlowWorkspace& workspace, const CaseKnownFlowLimits& limits,
                         core::ResourcePlan& memory, const std::stop_token& stop,
                         core::FailureDetails details) noexcept {
    return guarded<CaseKnownFlowJet>(evaluate, details, [&] {
        return known_flow_detail::Builder::run(requested, current, query, projection, workspace,
                                               limits, memory, stop, details);
    });
}
core::Result<CaseKnownFlowSnapshot>
snapshot_case_known_flow(const CaseKnownFlowReady& requested, a::UniformFlowWorkspace& workspace,
                         const CaseKnownFlowLimits& limits, core::ResourcePlan& memory,
                         const std::stop_token& stop, core::FailureDetails details) noexcept {
    return guarded<CaseKnownFlowSnapshot>(prepare, details, [&] {
        return known_flow_detail::Builder::save(requested, workspace, limits, memory, stop,
                                                details);
    });
}
core::Result<CaseKnownFlowReady>
restore_case_known_flow(const CaseKnownFlowSnapshot& saved, const CaseKnownFlowReady& current,
                        a::UniformFlowWorkspace& workspace, const CaseKnownFlowLimits& limits,
                        core::ResourcePlan& memory, const std::stop_token& stop,
                        core::FailureDetails details) noexcept {
    return guarded<CaseKnownFlowReady>(prepare, details, [&] {
        return known_flow_detail::Builder::restore(saved, current, workspace, limits, memory, stop,
                                                   details);
    });
}

} // namespace uvlas::coupling::detail
