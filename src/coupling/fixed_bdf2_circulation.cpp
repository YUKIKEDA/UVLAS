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
bool same_identity(const io::SchemaIdentity& x, const io::SchemaIdentity& y) {
    return x.source == y.source && x.contract_revision == y.contract_revision &&
           x.representation_revision == y.representation_revision;
}

template <class Sample>
bool same_case(const Sample& a, const Sample& b) {
    const auto x = a.source();
    const auto y = b.source();
    const auto ax = a.address();
    const auto ay = b.address();
    return same_identity(x, y) && ax.stamp == ay.stamp && ax.interval == ay.interval &&
           a.restart() == b.restart();
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

void verify_history(const CaseClockedFlowSample& clock, const DiscreteCirculation& result,
                    m::Work& work) {
    verify_ambient(clock, result.flow(), work);
}

void verify_known_time(const CaseClockedKnownFlowSample& clock, p::PchipInputStamp stamp,
                       p::RequestedSide side, p::RationalView time, m::Work& work) {
    const auto& ambient = clock.ambient();
    if (read(ambient.stamp()) != stamp || ambient.requested_side() != side ||
        !same_rational(read(clock.time()), time, work))
        reject(StaleGeneration);
}

void verify_history(const CaseClockedKnownFlowSample& clock, const DiscreteKnownCirculation& result,
                    m::Work& work) {
    if (!same_identity(clock.source(), read(result.case_identity())))
        reject(StaleGeneration);
    // The known solve samples every collocation with a single common flow query.
    const auto& sample = read(result.flow(0)).get();
    verify_known_time(clock, read(sample.stamp()), sample.requested_side(), read(sample.time()),
                      work);
}

a::PreparedTopologyOwner topology(const DiscreteCirculation& result) {
    return result.source();
}

a::PreparedTopologyOwner topology(const DiscreteKnownCirculation& result) {
    return read(result.source());
}

const numerics::linear::Solution<double>& retained(const DiscreteCirculation& result) {
    return result.solution();
}

const numerics::linear::Solution<double>& retained(const DiscreteKnownCirculation& result) {
    return read(result.solution()).get();
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

    using Destination = std::variant<CaseClockedFlowSample, CaseClockedKnownFlowSample>;
    std::optional<Destination> destination;
    std::optional<core::AcceptedStateId> parent;
    std::optional<a::PreparedTopologyOwner> source;
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
    return prepare_history(std::move(destination), previous, older, previous_solution,
                           older_solution, std::nullopt, maximum_integer_bits, maximum_steps,
                           memory, stop);
}

core::Result<FixedBdf2Circulation> FixedBdf2Circulation::prepare(
    CaseClockedKnownFlowSample&& destination, const CaseClockedKnownFlowSample& previous,
    const CaseClockedKnownFlowSample& older, const DiscreteKnownCirculation& previous_solution,
    const DiscreteKnownCirculation& older_solution, core::AcceptedStateId parent,
    std::uint64_t maximum_integer_bits, std::uint64_t maximum_steps, core::ResourcePlan& memory,
    const std::stop_token& stop) noexcept {
    return prepare_history(std::move(destination), previous, older, previous_solution,
                           older_solution, parent, maximum_integer_bits, maximum_steps, memory,
                           stop);
}

template <class Sample, class Solve>
core::Result<FixedBdf2Circulation> FixedBdf2Circulation::prepare_history(
    Sample&& destination, const Sample& previous, const Sample& older,
    const Solve& previous_solution, const Solve& older_solution,
    std::optional<core::AcceptedStateId> parent, std::uint64_t maximum_integer_bits,
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
                verify_history(previous, previous_solution, work);
                verify_history(older, older_solution, work);
                auto source = topology(previous_solution);
                if (!source.same_owner(topology(older_solution)))
                    reject(StaleGeneration);
                const auto twice_h =
                    m::scale(case_components::rational(h, maximum_integer_bits, work), 1);
                if (m::compare(twice_h.lower, {}) <= 0)
                    reject(InvalidInput);
                const auto& last = retained(previous_solution).values();
                const auto& before = retained(older_solution).values();
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
                owner.source.emplace(std::move(source));
                owner.parent = parent;
                owner.twice_step = twice_h;
                owner.destination.emplace(std::in_place_type<Sample>, std::move(destination));
                return FixedBdf2Circulation{std::move(data)};
            } catch (const core::Failure& error) {
                return std::unexpected(error);
            } catch (const m::ArithmeticFailure& error) {
                return std::unexpected(core::Failure{error.code, operation});
            }
        });
}

const FixedBdf2Circulation::Data&
FixedBdf2Circulation::validate_source(const a::PreparedTopologyOwner& current,
                                      core::NumericalGeneration generation, std::size_t size,
                                      m::Work& work) const {
    work.step();
    const auto storage = data_.values();
    if (storage.empty())
        reject(InvalidInput);
    const auto& owner = storage.front();
    if (!owner.source || !owner.destination)
        reject(InvalidInput);
    if (!current.same_owner(*owner.source) || owner.rows.size() != size ||
        generation == owner.history[0] || generation == owner.history[1])
        reject(StaleGeneration);
    return owner;
}

void FixedBdf2Circulation::bind(const a::PreparedTopologyOwner& current,
                                const ProjectedCaseUniformFlow& flow,
                                core::NumericalGeneration generation, std::size_t size,
                                m::Work& work) const {
    const auto& owner = validate_source(current, generation, size, work);
    const auto* clock = std::get_if<CaseClockedFlowSample>(&*owner.destination);
    if (clock == nullptr)
        reject(StaleGeneration);
    verify_ambient(*clock, flow, work);
}

void FixedBdf2Circulation::bind(const a::PreparedTopologyOwner& current,
                                const DiscreteKnownFlowInput& flow, core::AcceptedStateId parent,
                                core::NumericalGeneration generation, std::size_t size,
                                m::Work& work) const {
    const auto& owner = validate_source(current, generation, size, work);
    const auto* clock = std::get_if<CaseClockedKnownFlowSample>(&*owner.destination);
    if (clock == nullptr || owner.parent != parent ||
        !same_identity(clock->source(), read(flow.source.case_identity())) ||
        !same_identity(clock->source(), read(flow.current.case_identity())))
        reject(StaleGeneration);
    verify_known_time(*clock, flow.time.stamp, flow.time.side, flow.time.time, work);
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
