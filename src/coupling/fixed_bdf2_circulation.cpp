#include "fixed_bdf2_circulation.hpp"

#include "aero/p1_interval_private.hpp"
#include "case_component_precision_private.hpp"

#include <algorithm>
#include <cmath>

namespace uvlas::coupling::detail {
namespace a = aero::detail;
namespace m = a::p1;
namespace p = a::pchip;

namespace {
using enum core::FailureCode;
constexpr auto operation = core::Operation::EvaluateAero;

[[noreturn]] void reject(core::FailureCode reason) {
    throw m::ArithmeticFailure{reason};
}
template <class T>
T read(core::Result<T> result) {
    if (!result)
        throw result.error();
    return std::move(*result);
}

bool compare_words(std::span<const std::uint32_t> a, std::span<const std::uint32_t> b,
                   m::Work& work) {
    return std::ranges::equal(a, b, [&](auto lhs, auto rhs) {
        work.step();
        return lhs == rhs;
    });
}
bool same_rational(p::RationalView a, p::RationalView b, m::Work& work) {
    work.step();
    return a.sign == b.sign && compare_words(a.numerator, b.numerator, work) &&
           compare_words(a.denominator, b.denominator, work);
}
bool same_case(const CaseClockedFlowSample& a, const CaseClockedFlowSample& b) {
    const auto x = a.source();
    const auto y = b.source();
    const auto ax = a.address();
    const auto ay = b.address();
    return x.source == y.source && x.contract_revision == y.contract_revision &&
           x.representation_revision == y.representation_revision && ax.stamp == ay.stamp &&
           ax.interval == ay.interval && a.restart() == b.restart();
}

bool same_reference(const CirculationPhysicalSource& a, const CirculationPhysicalSource& b) {
    const auto& x = a.geometry;
    const auto& y = b.geometry;
    return a.topology.same_owner(b.topology) && x.fingerprint == y.fingerprint &&
           x.model == y.model && x.reference == y.reference && x.mapping == y.mapping &&
           x.connectivity == y.connectivity && x.policy == y.policy && x.kernel == y.kernel &&
           x.radius == y.radius;
}

void verify_ambient(const CaseClockedFlowSample& clock, const ProjectedCaseUniformFlow& flow,
                    m::Work& work) {
    const auto& x = clock.ambient().source();
    const auto& y = flow.source();
    if (x.requested_side != y.requested_side || read(x.jet.stamp()) != read(y.jet.stamp()) ||
        !same_rational(read(clock.time()), read(y.jet.time()), work))
        reject(StaleGeneration);
    for (std::size_t k = 0; k < 3; ++k)
        if (!same_rational(read(x.jet.velocity(k)), read(y.jet.velocity(k)), work) ||
            !same_rational(read(x.jet.time_partial(k)), read(y.jet.time_partial(k)), work))
            reject(StaleGeneration);
}

a::ScalarWords coefficient(const numerics::PrecisionBuffer<double>& words, std::size_t i) {
    if (const auto* binary = std::get_if<numerics::Binary64<double>>(&words.storage()))
        return a::SingleWord{binary->values.values()[i]};
    if (const auto* paired = std::get_if<numerics::TwoWord<double>>(&words.storage()))
        return a::DoubleWord{paired->hi.values()[i], paired->lo.values()[i]};
    reject(PrecisionUnsupported);
}

m::Interval lift(a::ScalarWords value, m::Work& work) {
    work.step();
    if (const auto* binary = std::get_if<a::SingleWord>(&value))
        return a::decode_p1_word(binary->value);
    const auto& paired = std::get<a::DoubleWord>(value);
    return m::add(a::decode_p1_word(paired.high), a::decode_p1_word(paired.low), work);
}
} // namespace

struct FixedBdf2Circulation::Data {
    struct Row {
        a::PhysicalFaceId face{0};
        a::ScalarWords previous{a::SingleWord{0}}, older{a::SingleWord{0}};
    };

    std::optional<CaseClockedFlowSample> destination;
    std::optional<CirculationPhysicalSource> source;
    std::array<core::NumericalGeneration, 2> history{core::NumericalGeneration{0},
                                                     core::NumericalGeneration{0}};
    core::PlannedBuffer<Row> rows;
    m::Interval twice_step{};
};

FixedBdf2Circulation::FixedBdf2Circulation(core::PlannedBuffer<Data>&& data) noexcept
    : data_(std::move(data)) {}
FixedBdf2Circulation::FixedBdf2Circulation(FixedBdf2Circulation&&) noexcept = default;
FixedBdf2Circulation& FixedBdf2Circulation::operator=(FixedBdf2Circulation&&) noexcept = default;
FixedBdf2Circulation::~FixedBdf2Circulation() = default;

core::Result<FixedBdf2Circulation> FixedBdf2Circulation::prepare(
    CaseClockedFlowSample&& destination, const CaseClockedFlowSample& previous,
    const CaseClockedFlowSample& older, const DiscreteCirculation& previous_solution,
    const DiscreteCirculation& older_solution, std::uint64_t maximum_integer_bits,
    std::uint64_t maximum_steps, core::ResourcePlan& memory, const std::stop_token& stop) noexcept {
    return core::guard_entry(
        operation, PrecisionUnresolved, {}, [&]() -> core::Result<FixedBdf2Circulation> {
            try {
                m::Work work{0, maximum_steps, stop};
                work.step();
                const auto node = destination.address().step;
                const auto h = read(destination.step_seconds());
                if (node < 2 || previous.address().step != node - 1 ||
                    older.address().step != node - 2 || !same_case(destination, previous) ||
                    !same_case(destination, older) ||
                    !same_rational(h, read(previous.step_seconds()), work) ||
                    !same_rational(h, read(older.step_seconds()), work))
                    reject(StaleGeneration);
                (void)read(destination.time());
                verify_ambient(previous, previous_solution.flow(), work);
                verify_ambient(older, older_solution.flow(), work);
                if (!same_reference(previous_solution.source(), older_solution.source()))
                    reject(StaleGeneration);
                const auto twice_h =
                    m::scale(case_components::rational(h, maximum_integer_bits, work), 1);
                if (m::compare(twice_h.lower, {}) <= 0)
                    reject(InvalidInput);
                const auto& last = previous_solution.solution().values();
                const auto& before = older_solution.solution().values();
                const auto size = last.shape();
                if (size.rows != before.shape().rows || before.shape().columns != 1 ||
                    size.rows <= 0 || size.columns != 1 || last.generation() == before.generation())
                    reject(StaleGeneration);

                auto data = read(memory.allocate<Data>(1, operation));
                const auto records = data.values();
                auto& owner = records.front();
                owner.rows = read(memory.allocate<Data::Row>(size.rows, operation));
                for (std::size_t i = 0; i < owner.rows.size(); ++i) {
                    work.step();
                    const auto pface = read(previous_solution.row(i)).face;
                    if (pface != read(older_solution.row(i)).face)
                        reject(StaleGeneration);
                    owner.rows.values()[i] = {pface, coefficient(last, i), coefficient(before, i)};
                }
                work.step();
                owner.history = {last.generation(), before.generation()};
                owner.source.emplace(previous_solution.source());
                owner.twice_step = twice_h;
                owner.destination.emplace(std::move(destination));
                return FixedBdf2Circulation{std::move(data)};
            } catch (const core::Failure& error) {
                return std::unexpected(error);
            } catch (const m::ArithmeticFailure& error) {
                return std::unexpected(core::Failure{error.code, operation});
            }
        });
}

void FixedBdf2Circulation::bind(const CirculationPhysicalSource& current,
                                const ProjectedCaseUniformFlow& flow,
                                core::NumericalGeneration generation, std::size_t size,
                                m::Work& work) const {
    work.step();
    const auto storage = data_.values();
    if (storage.empty())
        reject(InvalidInput);
    const auto& owner = storage.front();
    if (!owner.source || !owner.destination)
        reject(InvalidInput);
    if (!same_reference(current, *owner.source) || owner.rows.size() != size ||
        generation == owner.history[0] || generation == owner.history[1])
        reject(StaleGeneration);
    verify_ambient(*owner.destination, flow, work);
}

m::Interval FixedBdf2Circulation::value(std::size_t i, a::PhysicalFaceId face,
                                        const m::Interval& current, m::Work& work) const {
    work.step();
    const auto records = data_.values();
    const auto& owner = records.front();
    const auto rows = owner.rows.values();
    const auto& history = rows[i];
    if (history.face != face)
        reject(StaleGeneration);
    // Difference form cancels equal imported values before multiplication.
    const auto previous = lift(history.previous, work);
    const auto first = m::subtract(current, previous, work);
    const auto second = m::subtract(previous, lift(history.older, work), work);
    return m::divide(m::subtract(m::multiply(m::integer(3), first, work), second, work),
                     owner.twice_step, work);
}

m::Interval FixedBdf2Circulation::current_action(const m::Interval& seed, m::Work& work) const {
    return m::divide(m::multiply(m::integer(3), seed, work), data_.values().front().twice_step,
                     work);
}

} // namespace uvlas::coupling::detail
