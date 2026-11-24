#include "case_pg_support_motion.hpp"

#include "case_component_precision_private.hpp"
#include "case_pg_support_input.hpp"
#include "io/si_case_detail.hpp"
#include "structure/pg_reference_arithmetic.hpp"

#include <algorithm>
#include <cmath>
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

using RigidJet = std::array<RigidCoefficients, 3>;

void fill_jet(RigidJet& out, const CaseMotionEvaluation& motion, std::size_t orders,
              m::Work& work) {
    for (std::size_t d = 0; d < orders; ++d)
        out[d] = coefficients(motion, static_cast<p::JetOrder>(d), work);
}

void fill_rhs(const structure::PgRigidSupportRows& support, std::span<const RigidJet> jets,
              const CasePgSupportProjection& projection, const CasePgSupportLimits& limits,
              std::size_t orders, std::span<io::SiEnclosure> output, m::Work& work) {
    auto scratch = a::LocalRootScratch{};
    const auto map = support.rows();
    const auto original = support.chart().original().original_rows();
    for (std::size_t row = 0; row < map.size(); ++row) {
        work.step();
        const auto& source = map[row];
        std::array<m::Interval, axes> x{};
        for (std::size_t k = 0; k < axes; ++k)
            x[k] =
                rational(source.group ? source.reference_position[k] : original[row].prescribed[k],
                         limits.maximum_row_width_bits, work);
        const auto b = source.group
                           ? rational(source.translation, limits.maximum_row_width_bits, work)
                           : m::integer(0);
        for (std::size_t d = 0; d < orders; ++d)
            for (std::size_t axis = 0; axis < axes; ++axis) {
                auto value = d == 0 ? x[axis] : m::integer(0);
                if (source.group) {
                    const auto& jet = jets[*source.group][d];
                    value = m::multiply(jet.translation[axis], b, work);
                    for (std::size_t k = 0; k < axes; ++k)
                        value = m::add(value, m::multiply(jet.rotation[axis][k], x[k], work), work);
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
                if (support.group_count() > 1)
                    case_components::reject(InvalidInput);
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
                RigidJet jet{};
                fill_jet(jet, sample.motion, orders, work);
                fill_rhs(support, std::span{&jet, 1U}, projection, limits, orders, rhs.values(),
                         work);
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

namespace structural_support_detail {
namespace {
namespace s = io::schema;
using namespace structure;
constexpr auto preparation = core::Operation::PrepareModel;

template <class T, class F>
core::Result<T> guarded(core::Operation op, core::FailureDetails details, F f) noexcept {
    return core::guard_entry(op, PrecisionUnresolved, details, [&]() -> core::Result<T> {
        try {
            return f();
        } catch (const core::Failure& failure) {
            return std::unexpected(failure);
        } catch (const m::ArithmeticFailure& failure) {
            return std::unexpected(core::Failure{failure.code, op, details});
        }
    });
}

template <class T>
bool equal_elements(std::span<const T> a, std::span<const T> b, m::Work& work) {
    if (a.size() != b.size())
        return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        work.step();
        if (a[i] != b[i])
            return false;
    }
    return true;
}

bool equal_words(p::RationalView a, p::RationalView b, m::Work& work) {
    work.step();
    return a.sign == b.sign && equal_elements(a.numerator, b.numerator, work) &&
           equal_elements(a.denominator, b.denominator, work);
}

bool equal_array(std::span<const p::RationalView> a, std::span<const p::RationalView> b,
                 m::Work& work) {
    if (a.size() != b.size())
        return false;
    for (std::size_t i = 0; i < a.size(); ++i)
        if (!equal_words(a[i], b[i], work))
            return false;
    return true;
}

bool equal_recipe(const CaseMotionSeriesSnapshotView& a, const CaseMotionSeriesSnapshotView& b,
                  m::Work& work) {
    if (a.index() != b.index())
        return false;
    if (const auto* f = std::get_if<CaseMotionFormulaSnapshotView>(&a)) {
        const auto& g = std::get<CaseMotionFormulaSnapshotView>(b);
        if (f->domain.has_value() != g.domain.has_value())
            return false;
        if (f->domain && g.domain && !equal_array(*f->domain, *g.domain, work))
            return false;
        return equal_array(f->coefficients, g.coefficients, work);
    }
    const auto& f = std::get<p::PchipSnapshotView>(a);
    const auto& g = std::get<p::PchipSnapshotView>(b);
    const auto& x = f.input;
    const auto& y = g.input;
    if (f.formula_version != g.formula_version || x.stamp != y.stamp ||
        x.samples.components != y.samples.components ||
        x.boundary.policy.index() != y.boundary.policy.index())
        return false;
    if (const auto* period = std::get_if<p::PeriodicBoundary>(&x.boundary.policy))
        if (!equal_words(period->period, std::get<p::PeriodicBoundary>(y.boundary.policy).period,
                         work))
            return false;
    return std::ranges::equal(x.dimensions, y.dimensions,
                              [&](const auto& l, const auto& r) {
                                  work.step();
                                  return l.powers == r.powers;
                              }) &&
           equal_words(x.boundary.begin, y.boundary.begin, work) &&
           equal_words(x.boundary.end, y.boundary.end, work) &&
           equal_array(x.events, y.events, work) &&
           equal_array(x.samples.times, y.samples.times, work) &&
           equal_array(x.samples.values, y.samples.values, work) &&
           equal_array(f.coefficients, g.coefficients, work);
}

bool equal_motion(const CaseMotionSnapshotView& a, const CaseMotionSnapshotView& b, m::Work& work) {
    return a.formula_version == b.formula_version &&
           io::detail::same_identity(a.frame.source, b.frame.source) &&
           io::detail::same_identity(a.time.source, b.time.source) &&
           a.time.stamp == b.time.stamp && a.time.pi_power == b.time.pi_power &&
           equal_array(a.frame.components, b.frame.components, work) &&
           equal_words(a.time.begin, b.time.begin, work) &&
           equal_words(a.time.end, b.time.end, work) &&
           equal_array(a.time.events, b.time.events, work) &&
           equal_recipe(a.time.translation, b.time.translation, work) &&
           equal_recipe(a.time.rotation, b.time.rotation, work);
}

s::SupportTraceConstraint from_case(const s::Model& model, PgConstraintId id, m::Work& work) {
    for (std::uint64_t i = 0; i < model.supports.size(); ++i) {
        work.step();
        const auto declared = take(model.supports.at(static_cast<core::Index>(i))).record();
        if (const auto* motion = std::get_if<s::SupportTraceConstraint>(&declared))
            if (motion->id.unsigned_value() == id.value())
                return *motion;
    }
    case_components::reject(InvalidInput);
}

void require_whole_edge(const s::SupportTraceConstraint& support, const s::Model& model,
                        const PgConstraintSource& source, m::Work& work) {
    const auto region = support.region.record();
    const auto target = region.target.record();
    if (target.kind.selected() != decltype(target.kind.selected())::ItemShell ||
        target.id.unsigned_value() != source.left.value())
        case_components::reject(StaleGeneration);
    if (region.kind.selected() != decltype(region.kind.selected())::ItemEdge ||
        region.side.selected() != decltype(region.side.selected())::ItemInterior ||
        support.components.size() != 3)
        case_components::reject(Unsupported);
    std::array<bool, 3> seen{};
    for (core::Index i = 0; i < 3; ++i) {
        work.step();
        const auto axis = take(support.components.at(i)).unsigned_value();
        if (axis >= seen.size() || seen[axis])
            case_components::reject(InvalidInput);
        seen[axis] = true;
    }

    for (std::uint64_t i = 0; i < model.shells.size(); ++i) {
        work.step();
        const auto shell = take(model.shells.at(static_cast<core::Index>(i))).record();
        if (shell.id.unsigned_value() == source.left.value()) {
            if (shell.geometry.unsigned_value() != region.geometry.unsigned_value())
                case_components::reject(StaleGeneration);
            return;
        }
    }
    case_components::reject(StaleGeneration);
}

core::Result<void> validate_region(const PgGlobalConstraints& original,
                                   const PgConstraintSource& source, const PgPreparedPatch& patch,
                                   std::uint64_t boundary, const PgTraceLimits& trace_limits,
                                   core::Workspace& workspace, core::ResourcePlan& plan,
                                   bool stationary = false) noexcept {
    using namespace structure::detail;
    return partition_entry([&]() -> core::Result<void> {
        ReferenceArithmetic math{trace_limits.arithmetic, plan};
        math.cancel();
        math.admit();
        const auto basis = partition_take(patch.conditioning().reference_basis());
        const auto& partition = partition_take(basis->source()).partition;
        const auto traces = partition_take(partition.edge_traces());
        const auto descriptors = partition_take(partition.sources());
        const auto overlap_count = static_cast<std::size_t>(source.overlaps);
        auto used = partition_take(
            plan.allocate<bool>(partition_index(overlap_count), partition_operation));
        const auto intersections = original.overlaps().subspan(
            static_cast<std::size_t>(source.first_overlap), overlap_count);
        std::size_t expected_count = 0;
        for (const auto& trace : traces) {
            math.visit(1);
            const auto& descriptor = descriptors[static_cast<std::size_t>(trace.source)].descriptor;
            const auto* named = std::get_if<PgPartitionNamedSource>(&descriptor);
            if (named == nullptr || !std::cmp_equal(named->boundary, boundary))
                continue;
            ++expected_count;
            const auto match = std::ranges::find_if(intersections, [&](const auto& overlap) {
                math.visit(1);
                return overlap.left_edge == trace.edge;
            });
            if (match == intersections.end() || match->left_layer)
                partition_fail(StaleGeneration);
            const auto i = static_cast<std::size_t>(match - intersections.begin());
            if (used.values()[i])
                partition_fail(StaleGeneration);
            used.values()[i] = true;
        }
        if (expected_count == 0 || expected_count != overlap_count)
            partition_fail(StaleGeneration);

        // Re-form one full named segment at a time in the existing row order.
        const auto zero = ExactRational{};
        const auto one = math.constant("1");
        const auto identity = partition_take(patch.info()).identity;
        const auto layout = original.layout();
        const auto patch_row = std::ranges::find_if(layout, [&](const auto& x) {
            math.visit(1);
            return x.id == source.left;
        });
        if (patch_row == layout.end())
            partition_fail(StaleGeneration);
        std::size_t matched_rows = 0;
        const auto row_count = static_cast<std::size_t>(source.rows);
        for (std::size_t i = 0; i < overlap_count; ++i) {
            const auto& overlap = intersections[i];
            const auto it = std::ranges::find_if(traces, [&](const auto& edge) {
                math.visit(1);
                const auto& descriptor =
                    descriptors[static_cast<std::size_t>(edge.source)].descriptor;
                const auto* named = std::get_if<PgPartitionNamedSource>(&descriptor);
                return named != nullptr && std::cmp_equal(named->boundary, boundary) &&
                       edge.edge == overlap.left_edge;
            });
            if (it == traces.end())
                partition_fail(StaleGeneration);
            const std::array spans{PgTraceEdgeSpan{overlap.left_edge,
                                                   overlap.first,
                                                   overlap.last,
                                                   it->canonical_direction ? zero : one,
                                                   it->canonical_direction ? one : zero,
                                                   {}}};
            const auto expected = partition_take(
                PgFullTrace::prepare_support(patch, spans, identity, patch.reference_coordinates(),
                                             trace_limits, workspace, plan));
            for (const auto& row : expected.original_rows()) {
                math.visit(1);
                if (matched_rows == row_count)
                    partition_fail(StaleGeneration);
                const auto index =
                    partition_index(static_cast<std::uint64_t>(source.first_row) + matched_rows);
                ++matched_rows;
                const auto actual = partition_take(original.row(index));
                const auto& recorded = original.original_rows()[static_cast<std::size_t>(index)];
                if (recorded.field != row.field || recorded.power != row.power ||
                    recorded.source_overlap != partition_index(i))
                    partition_fail(StaleGeneration);
                if (stationary)
                    for (std::size_t axis = 0; axis < axes; ++axis)
                        if (math.compare(recorded.prescribed[axis], row.prescribed[axis]) != 0)
                            partition_fail(StaleGeneration);
                if (row.unique_row < 0) {
                    if (!actual.empty())
                        partition_fail(StaleGeneration);
                    continue;
                }
                const auto base = partition_take(expected.unique_row(row.unique_row));
                if (actual.size() != base.size())
                    partition_fail(StaleGeneration);
                for (std::size_t j = 0; j < base.size(); ++j) {
                    math.visit(1);
                    const auto coefficient = math.mul(row.scale, base[j].coefficient);
                    if (actual[j].column != patch_row->first + base[j].column ||
                        math.compare(actual[j].coefficient, coefficient) != 0)
                        partition_fail(StaleGeneration);
                }
            }
        }
        if (matched_rows != row_count)
            partition_fail(StaleGeneration);

        return {};
    });
}
} // namespace

struct Builder {
    static CaseStructuralSupport
    prepare(const io::SiCase& si, PgSparseChart&& chart,
            std::span<const PgConstraintPatchBinding> patches, CasePgSupportSelection selection,
            p::PchipInputStamp stamp, core::Workspace& pg_workspace,
            p::CurveWorkspace& motion_workspace, const CaseStructuralSupportLimits& limits,
            core::ResourcePlan& plan, const std::stop_token& stop, core::FailureDetails details) {
        m::Work work{0, limits.maximum_group_work, stop};
        work.step();
        const auto ids = selection.moving;
        const auto stationary = selection.stationary;
        if (ids.empty())
            case_components::reject(InvalidInput);
        if (ids.size() > limits.maximum_group_supports ||
            stationary.size() > limits.maximum_group_supports - ids.size() ||
            !std::in_range<core::Index>(ids.size()) ||
            !std::in_range<core::Index>(stationary.size()))
            case_components::reject(ResourceExceeded);
        const auto& timing = limits.evaluation.clock;
        auto grid = take(prepare_case_reference_clock(
            si, stamp, timing.flow.clock_and_wake, timing.maximum_event_work, plan, stop, details));
        const auto intervals = grid.intervals.values();
        if (intervals.empty())
            case_components::reject(InvalidInput);
        const auto root = take(take(si.source_case()).get().root());
        const auto model = root.model.record();
        auto order =
            take(plan.allocate<std::uint64_t>(static_cast<core::Index>(ids.size()), preparation));
        for (std::size_t i = 0; i < ids.size(); ++i) {
            work.step();
            order.values()[i] = ids[i].value();
        }
        std::ranges::sort(order.values(), [&](auto a, auto b) {
            work.step();
            return a < b;
        });
        const auto capacity = std::min(ids.size(), limits.maximum_motion_groups);
        if (capacity == 0)
            case_components::reject(ResourceExceeded);
        auto motions = take(plan.allocate<std::optional<CaseMotionReady>>(
            static_cast<core::Index>(capacity), preparation));
        auto archives = take(plan.allocate<std::optional<CaseMotionSnapshot>>(
            static_cast<core::Index>(capacity), preparation));
        auto bindings = take(plan.allocate<PgRigidSupportBinding>(
            static_cast<core::Index>(ids.size()), preparation));
        std::size_t group_count = 0;
        const auto& original = chart.original();
        const auto sources = original.sources();
        const auto in_group = [&](std::uint64_t id) {
            return std::ranges::binary_search(order.values(), id, [&](auto a, auto b) {
                work.step();
                return a < b;
            });
        };
        const auto is_stationary = [&](std::uint64_t id) {
            return std::ranges::any_of(stationary, [&](auto entry) {
                work.step();
                return entry.value() == id;
            });
        };
        for (const auto& source : sources) {
            work.step();
            if (source.kind == PgTraceKind::Support && !in_group(source.id.value()) &&
                !is_stationary(source.id.value()))
                case_components::reject(Unsupported);
        }
        for (std::uint64_t i = 0; i < model.supports.size(); ++i) {
            work.step();
            const auto declared = take(model.supports.at(static_cast<core::Index>(i))).record();
            if (const auto* support = std::get_if<s::SupportTraceConstraint>(&declared)) {
                const auto target = support->region.record().target.record();
                if (target.kind.selected() != decltype(target.kind.selected())::ItemShell)
                    continue;
                const auto shell = target.id.unsigned_value();
                const bool in_layout = std::ranges::any_of(patches, [&](const auto& patch) {
                    work.step();
                    return patch.id.value() == shell;
                });
                if (in_layout && !in_group(support->id.unsigned_value()) &&
                    !is_stationary(support->id.unsigned_value()))
                    case_components::reject(Unsupported);
            }
        }
        auto trace_limits = limits.region_trace;
        trace_limits.arithmetic.cancellation = stop;
        const auto validate_support = [&](PgConstraintId id, bool fixed) {
            const auto declaration = from_case(model, id, work);
            const auto source = std::ranges::find_if(sources, [&](const auto& x) {
                work.step();
                return x.kind == PgTraceKind::Support && x.id == id;
            });
            if (source == sources.end())
                case_components::reject(StaleGeneration);
            require_whole_edge(declaration, model, *source, work);
            const auto patch = std::ranges::find_if(patches, [&](const auto& x) {
                work.step();
                return x.id == source->left;
            });
            if (patch == patches.end())
                case_components::reject(StaleGeneration);
            const auto checked =
                validate_region(original, *source, patch->patch,
                                declaration.region.record().boundary.unsigned_value(), trace_limits,
                                pg_workspace, plan, fixed);
            if (!checked)
                throw checked.error();
            return declaration;
        };
        for (std::size_t i = 0; i < stationary.size(); ++i) {
            work.step();
            const auto id = stationary[i];
            if (in_group(id.value()) ||
                std::ranges::any_of(stationary.first(i), [&](auto previous) {
                    work.step();
                    return previous == id;
                }))
                case_components::reject(InvalidInput);
            const auto declaration = validate_support(id, true);
            // The event count and at most six component lookups.
            for (std::size_t lookup = 0; lookup < 7; ++lookup)
                work.step();
            if (!take(stationary_case_pg_motion(si, declaration.motion.record())))
                case_components::reject(Unsupported);
            const auto stationary_frame = take(CaseMotionFrame::prepare(
                si, declaration.motion, timing.motion.frame, plan, stop, details));
            static_cast<void>(stationary_frame);
        }
        for (std::size_t i = 0; i < ids.size(); ++i) {
            work.step();
            const auto id = PgConstraintId{order.values()[i]};
            if (i != 0 && id.value() == order.values()[i - 1])
                case_components::reject(InvalidInput);
            const auto declaration = validate_support(id, false);
            auto selected = take(prepare_case_motion(
                si, declaration.motion, {intervals.front().begin, intervals.back().end}, stamp,
                motion_workspace, timing.motion, plan, stop, details));
            auto snapshot = take(snapshot_case_motion(selected, motion_workspace, timing.motion,
                                                      plan, stop, details));
            const auto view = take(snapshot.view());
            const auto occupied = archives.values().first(group_count);
            const auto match = std::ranges::find_if(occupied, [&](const auto& known) {
                work.step();
                return equal_motion(take(known->view()), view, work);
            });
            const auto group = static_cast<std::size_t>(match - occupied.begin());
            if (match == occupied.end()) {
                if (group_count == capacity)
                    case_components::reject(ResourceExceeded);
                motions.values()[group].emplace(std::move(selected));
                archives.values()[group].emplace(std::move(snapshot));
                ++group_count;
            }
            bindings.values()[i] = {id, group};
        }
        auto exact = limits.exact;
        exact.arithmetic.cancellation = stop;
        auto support = take(PgRigidSupportRows::prepare_groups(
            std::move(chart), patches, bindings.values(), group_count, exact, pg_workspace, plan));
        return {std::move(grid), std::move(motions), std::move(support)};
    }

    static CaseStructuralSupportSample
    evaluate(const CaseStructuralSupport& field, const CaseStructuralSupport& current,
             CaseFlowClockQuery address, p::JetOrder order,
             const CasePgSupportProjection& projection, p::CurveWorkspace& motion_workspace,
             const CasePgSupportLimits& limits, core::ResourcePlan& plan,
             const std::stop_token& stop, core::FailureDetails details) {
        m::Work work{0, limits.maximum_row_work, stop};
        work.step();
        const auto orders = static_cast<std::size_t>(order) + 1;
        if (orders > 3 || projection.motion.size() != 12 * orders ||
            projection.probe.size() != axes * orders ||
            projection.original_rhs.size() != axes * orders)
            case_components::reject(InvalidInput);
        const auto id = take(field.support_.chart().identity());
        const auto map = field.support_.rows(), live = current.support_.rows();
        const auto groups = field.support_.group_count();
        if (id != take(current.support_.chart().identity()) || map.size() != live.size() ||
            current.support_.group_count() != groups || groups > field.motions_.size() ||
            groups > current.motions_.size())
            case_components::reject(StaleGeneration);
        for (std::size_t i = 0; i < map.size(); ++i) {
            work.step();
            if (map[i].group != live[i].group)
                case_components::reject(StaleGeneration);
        }
        if (map.size() > limits.maximum_original_rows ||
            map.size() > static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max()) /
                             (axes * orders))
            case_components::reject(ResourceExceeded);
        auto node = take(evaluate_case_reference_clock(field.clock_, current.clock_, address,
                                                       limits.clock.flow.clock_and_wake, plan, stop,
                                                       details));
        const auto probes = reference_probes(projection.probe);
        const auto group_size = static_cast<core::Index>(groups);
        auto result =
            take(plan.allocate<std::optional<CaseMotionEvaluation>>(group_size, operation));
        auto jets = take(plan.allocate<RigidJet>(group_size, operation));
        for (std::size_t group = 0; group < groups; ++group) {
            work.step();
            const auto& motion = field.motions_.values()[group];
            const auto& current_motion = current.motions_.values()[group];
            if (!motion || !current_motion)
                case_components::reject(StaleGeneration);
            const auto& published = result.values()[group].emplace(take(evaluate_case_motion(
                *motion, *current_motion, {address.stamp, node.seconds, order, address.side},
                {projection.representation, projection.motion, probes}, motion_workspace,
                limits.clock.motion, plan, stop, details)));
            fill_jet(jets.values()[group], published, orders, work);
        }
        auto rhs = take(plan.allocate<io::SiEnclosure>(
            static_cast<core::Index>(map.size() * axes * orders), operation));
        fill_rhs(field.support_, jets.values(), projection, limits, orders, rhs.values(), work);
        work.step();
        auto row_groups = take(plan.allocate<std::optional<std::size_t>>(
            static_cast<core::Index>(map.size()), operation));
        for (std::size_t i = 0; i < map.size(); ++i) {
            work.step();
            row_groups.values()[i] = map[i].group;
        }
        return {std::move(node), std::move(result),     id,         field.source(),
                std::move(rhs),  std::move(row_groups), map.size(), orders};
    }
};
} // namespace structural_support_detail

core::Result<std::reference_wrapper<const CaseMotionEvaluation>>
CaseStructuralSupportSample::motion(std::size_t group) const& noexcept {
    if (group >= motions_.size())
        return std::unexpected(core::Failure{InvalidInput, operation});
    const auto& result = motions_.values()[group];
    if (!result)
        return std::unexpected(core::Failure{StaleGeneration, operation});
    return std::cref(*result);
}

core::Result<CaseStructuralSupport> prepare_case_structural_pg_support(
    const io::SiCase& si, structure::PgSparseChart&& chart,
    std::span<const structure::PgConstraintPatchBinding> patches,
    std::span<const structure::PgConstraintId> ids, p::PchipInputStamp stamp,
    core::Workspace& pg_workspace, p::CurveWorkspace& motion_workspace,
    const CaseStructuralSupportLimits& limits, core::ResourcePlan& plan,
    const std::stop_token& stop, core::FailureDetails details) noexcept {
    return prepare_case_structural_pg_support(si, std::move(chart), patches,
                                              CasePgSupportSelection{ids, {}}, stamp, pg_workspace,
                                              motion_workspace, limits, plan, stop, details);
}

core::Result<CaseStructuralSupport> prepare_case_structural_pg_support(
    const io::SiCase& si, structure::PgSparseChart&& chart,
    std::span<const structure::PgConstraintPatchBinding> patches, CasePgSupportSelection selection,
    p::PchipInputStamp stamp, core::Workspace& pg_workspace, p::CurveWorkspace& motion_workspace,
    const CaseStructuralSupportLimits& limits, core::ResourcePlan& plan,
    const std::stop_token& stop, core::FailureDetails details) noexcept {
    using structural_support_detail::Builder;
    return structural_support_detail::guarded<CaseStructuralSupport>(
        core::Operation::PrepareModel, details, [&] {
            return Builder::prepare(si, std::move(chart), patches, selection, stamp, pg_workspace,
                                    motion_workspace, limits, plan, stop, details);
        });
}

core::Result<CaseStructuralSupportSample> evaluate_case_structural_pg_support(
    const CaseStructuralSupport& field, const CaseStructuralSupport& current,
    CaseFlowClockQuery address, p::JetOrder order, const CasePgSupportProjection& projection,
    p::CurveWorkspace& motion_workspace, const CasePgSupportLimits& limits,
    core::ResourcePlan& plan, const std::stop_token& stop, core::FailureDetails details) noexcept {
    using structural_support_detail::Builder;
    return structural_support_detail::guarded<CaseStructuralSupportSample>(operation, details, [&] {
        return Builder::evaluate(field, current, address, order, projection, motion_workspace,
                                 limits, plan, stop, details);
    });
}

namespace structural_support_detail {
namespace {
double positive_upper(const m::Bound& value) {
    if (value.mantissa < 0)
        case_components::reject(InvalidInput);
    if (value.mantissa == 0)
        return 0.0;
    const auto bits = static_cast<int>(boost::multiprecision::msb(value.mantissa)) + 1;
    const auto top = static_cast<std::int64_t>(value.power) + bits - 1;
    if (top > 1023)
        case_components::reject(PrecisionUnresolved);
    if (top < -1074)
        return std::numeric_limits<double>::denorm_min();
    const auto shift = static_cast<unsigned>(std::max({0, bits - 53, -1074 - value.power}));
    auto significand = m::Integer{value.mantissa >> shift};
    if ((significand << shift) != value.mantissa)
        ++significand;
    const auto word =
        std::scalbn(significand.convert_to<double>(), value.power + static_cast<int>(shift));
    if (!std::isfinite(word))
        case_components::reject(PrecisionUnresolved);
    return word;
}

double unscaled_word(double value, std::int64_t exponent) {
    const auto e = m::exponent(exponent);
    const auto scaled = std::scalbn(value, e);
    if (!std::isfinite(scaled) || m::compare(m::scale(a::decode_p1_word(value), e).lower,
                                             a::decode_p1_word(scaled).lower) != 0)
        case_components::reject(PrecisionUnresolved);
    return scaled;
}

m::Interval absolute_words(const numerics::PrecisionBuffer<double>& values, std::size_t i,
                           m::Work& work) {
    if (const auto* high = std::get_if<numerics::Binary64<double>>(&values.storage()))
        return m::absolute(a::decode_p1_word(high->values.values()[i]));
    const auto& words = std::get<numerics::TwoWord<double>>(values.storage());
    return m::absolute(m::add(a::decode_p1_word(words.hi.values()[i]),
                              a::decode_p1_word(words.lo.values()[i]), work));
}

double checked_upper(const m::Interval& bound, double allowance) {
    if (!std::isfinite(allowance))
        case_components::reject(Nonfinite);
    if (allowance < 0)
        case_components::reject(InvalidInput);
    if (m::compare(bound.upper, a::decode_p1_word(allowance).lower) > 0)
        case_components::reject(PrecisionUnresolved);
    return positive_upper(bound.upper);
}
} // namespace

struct Retention {
    static CaseStructuralRetainedRhs apply(const CaseStructuralSupportSample& sample,
                                           p::JetOrder order, core::NumericalGeneration generation,
                                           std::uint64_t maximum_work, core::ResourcePlan& plan,
                                           const std::stop_token& stop) {
        m::Work work{0, maximum_work, stop};
        work.step();
        if (!a::p1_environment_supported())
            case_components::reject(PrecisionUnsupported);
        const auto d = static_cast<std::size_t>(order);
        const auto rows = sample.rows();
        if (sample.orders() == 0 || sample.orders() > 3 || d >= sample.orders() ||
            rows > static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max()) /
                       (axes * sample.orders()) ||
            sample.original_rhs().size() != rows * axes * sample.orders())
            case_components::reject(InvalidInput);
        const auto count = static_cast<core::Index>(rows * axes);
        auto high = take(plan.allocate<double>(count, operation));
        auto low = take(plan.allocate<double>(count, operation));
        auto radii = take(plan.allocate<double>(count, operation));
        double max_radius = 0;
        bool paired = false;
        for (std::size_t row = 0; row < rows; ++row)
            for (std::size_t component = 0; component < axes; ++component) {
                work.step();
                const auto out = component * rows + row;
                const auto& in = sample.original_rhs()[(d * rows + row) * axes + component];
                std::visit(
                    [&](const auto& word) {
                        if constexpr (std::is_same_v<std::remove_cvref_t<decltype(word)>,
                                                     io::ScaledBinary64>) {
                            high.values()[out] = unscaled_word(word.center, word.exponent2);
                        } else {
                            high.values()[out] = unscaled_word(word.high, word.exponent2);
                            low.values()[out] = unscaled_word(word.low, word.exponent2);
                            paired = true;
                        }
                        radii.values()[out] = positive_upper(
                            m::scale(a::decode_p1_word(word.radius), m::exponent(word.exponent2))
                                .upper);
                        max_radius = std::max(max_radius, radii.values()[out]);
                    },
                    in);
            }
        const auto n = static_cast<core::Index>(rows);
        const numerics::DenseShape shape{n, 3, n};
        auto first = take(numerics::DenseReal::copy(shape, high.values(), plan, operation));
        using Buffer = numerics::PrecisionBuffer<double>;
        auto centre =
            paired ? take(Buffer::two_word(
                         std::move(first),
                         take(numerics::DenseReal::copy(shape, low.values(), plan, operation)),
                         generation))
                   : take(Buffer::binary64(std::move(first), generation));
        work.step();
        return {sample.chart(), std::move(centre), std::move(radii), max_radius};
    }
};
} // namespace structural_support_detail

core::Result<CaseStructuralRetainedRhs>
retain_case_structural_pg_rhs(const CaseStructuralSupportSample& sample, p::JetOrder order,
                              core::NumericalGeneration generation, std::uint64_t maximum_work,
                              core::ResourcePlan& plan, const std::stop_token& stop,
                              core::FailureDetails details) noexcept {
    using structural_support_detail::Retention;
    return structural_support_detail::guarded<CaseStructuralRetainedRhs>(operation, details, [&] {
        return Retention::apply(sample, order, generation, maximum_work, plan, stop);
    });
}

core::Result<structure::PgChartCandidateLift> lift_case_structural_pg_support(
    structure::PgChartCandidate& chart, structure::PgChartCandidateInput independent,
    const CaseStructuralRetainedRhs& rhs, const structure::PgChartCandidateSolve& request,
    core::Workspace& workspace, core::ResourcePlan& plan, core::FailureDetails details) noexcept {
    return structural_support_detail::guarded<structure::PgChartCandidateLift>(
        operation, details, [&] {
            m::Work work{0, 1, request.cancellation};
            auto result = take(
                chart.lift(independent, {rhs.chart(), rhs.centre()}, request, workspace, plan));
            const auto bound = m::add(a::decode_p1_word(result.original_residual_upper),
                                      a::decode_p1_word(rhs.maximum_radius()), work);
            result.original_residual_upper = structural_support_detail::checked_upper(
                bound, request.maximum_absolute_original_defect);
            return result;
        });
}

core::Result<double> work_case_structural_pg_support(
    const structure::PgChartCandidate& chart, structure::PgChartCandidateInput force,
    structure::PgChartCandidateInput independent, const CaseStructuralRetainedRhs& rhs,
    const structure::PgChartCandidateLift& lift, const structure::PgChartCandidateDual& dual,
    double tolerance, std::uint64_t max_work, const std::stop_token& stop,
    core::FailureDetails details) noexcept {
    return structural_support_detail::guarded<double>(operation, details, [&] {
        m::Work work{0, max_work, stop};
        work.step();
        const auto defect = take(chart.work_defect(force, independent, {rhs.chart(), rhs.centre()},
                                                   lift, dual, tolerance, max_work, stop));
        auto bound = a::decode_p1_word(defect);
        const auto shape = dual.original_multipliers.shape();
        const auto rows = static_cast<std::size_t>(shape.rows);
        if (rhs.radii().size() != rows * axes)
            case_components::reject(InvalidInput);
        for (std::size_t axis = 0; axis < axes; ++axis)
            for (std::size_t row = 0; row < rows; ++row) {
                work.step();
                const auto i = axis * rows + row;
                const auto coefficient = structural_support_detail::absolute_words(
                    dual.original_multipliers,
                    axis * static_cast<std::size_t>(shape.leading_dimension) + row, work);
                bound = m::add(
                    bound, m::multiply(coefficient, a::decode_p1_word(rhs.radii()[i]), work), work);
            }
        return structural_support_detail::checked_upper(bound, tolerance);
    });
}

core::Result<CaseStructuralExactLift> lift_case_structural_pg_support_exact(
    const CaseStructuralSupport& field, const CaseStructuralSupportSample& sample,
    p::JetOrder order, std::span<const structure::PgCartesianExactPosition> independent,
    const structure::PgRigidSupportLimits& limits, core::Workspace& workspace,
    core::ResourcePlan& plan) noexcept {
    using namespace structure;
    using namespace structure::detail;
    return partition_entry([&]() -> core::Result<CaseStructuralExactLift> {
        ReferenceArithmetic math{limits.arithmetic, plan};
        math.cancel();
        if (!a::p1_environment_supported())
            partition_fail(PrecisionUnsupported);
        const auto& chart = field.support().chart();
        const auto identity = partition_take(chart.identity());
        if (identity != sample.chart() ||
            !io::detail::same_identity(field.source(), sample.source()))
            partition_fail(StaleGeneration);
        const auto map = field.support().rows();
        const auto d = static_cast<std::size_t>(order);
        const auto source = sample.original_rhs();
        if (sample.rows() != map.size() || d >= sample.orders() || sample.orders() > 3 ||
            source.size() != partition_mul(partition_mul(map.size(), axes), sample.orders()))
            partition_fail(InvalidInput);
        const auto groups = field.support().group_count();
        const auto selected = sample.row_groups();
        if (selected.size() != map.size() || sample.motion_count() != groups)
            partition_fail(StaleGeneration);
        for (std::size_t i = 0; i < map.size(); ++i) {
            math.visit(1);
            if (selected[i] != map[i].group)
                partition_fail(StaleGeneration);
        }
        if (map.size() > limits.maximum_original_rows)
            partition_fail(ResourceExceeded);
        const auto growing = partition_mul(axes, map.size());
        const auto group_rationals = partition_mul(axes * (axes + 1), groups);
        math.admit(
            partition_add(reference_live_rationals, partition_add(growing, group_rationals)));

        const auto centre = [&](const io::SiEnclosure& input, bool include_radius) {
            return std::visit(
                [&](const auto& w) {
                    const auto exact = [&](double stored) {
                        math.charge();
                        auto word = partition_take(
                            ExactRational::from_binary64(stored, math.leaf_limits, plan));
                        if (word.sign() == 0 || w.exponent2 == 0)
                            return word;
                        const auto exponent = static_cast<std::int64_t>(m::exponent(w.exponent2));
                        auto bits = static_cast<std::uint64_t>(exponent < 0 ? -exponent : exponent);
                        auto factor = math.constant("2");
                        while (bits != 0) {
                            math.visit(1);
                            if (bits % 2 != 0)
                                word =
                                    exponent < 0 ? math.div(word, factor) : math.mul(word, factor);
                            bits /= 2;
                            if (bits != 0)
                                factor = math.square(factor);
                        }
                        return word;
                    };
                    if (include_radius)
                        return exact(w.radius);
                    if constexpr (std::is_same_v<std::remove_cvref_t<decltype(w)>,
                                                 io::ScaledBinary64>)
                        return exact(w.center);
                    else
                        return math.add(exact(w.high), exact(w.low));
                },
                input);
        };
        const auto compute = [&]() -> core::Result<CaseStructuralExactLift> {
            auto rhs = partition_take(plan.allocate<PgCartesianExactPosition>(
                partition_index(map.size()), partition_operation));
            ExactRational maximum_error;
            {
                const auto lease = partition_take(workspace.acquire());
                static_cast<void>(lease);
                struct GroupJet {
                    PgCartesianExactPosition translation;
                    std::array<PgCartesianExactPosition, axes> rotation;
                };
                auto jets = partition_take(
                    plan.allocate<GroupJet>(partition_index(groups), partition_operation));
                for (std::size_t group = 0; group < groups; ++group) {
                    math.visit(1);
                    const auto& motion = partition_take(sample.motion(group)).get();
                    auto& jet = jets.values()[group];
                    for (std::size_t axis = 0; axis < axes; ++axis) {
                        jet.translation[axis] =
                            centre(partition_take(motion.point(0, axis, order)), false);
                        for (std::size_t j = 0; j < axes; ++j)
                            jet.rotation[axis][j] = math.sub(
                                centre(partition_take(motion.point(j + 1, axis, order)), false),
                                jet.translation[axis]);
                    }
                }
                const auto original = chart.original().original_rows();
                for (std::size_t row = 0; row < map.size(); ++row) {
                    math.visit(1);
                    for (std::size_t axis = 0; axis < axes; ++axis) {
                        auto& output = rhs.values()[row][axis];
                        if (map[row].group) {
                            const auto& jet = jets.values()[*map[row].group];
                            output = math.mul(jet.translation[axis], map[row].translation);
                            for (std::size_t j = 0; j < axes; ++j)
                                output = math.add(output, math.mul(jet.rotation[axis][j],
                                                                   map[row].reference_position[j]));
                        } else if (d == 0) {
                            output = math.copy(original[row].prescribed[axis]);
                        }
                        const auto& approximated = source[(d * map.size() + row) * axes + axis];
                        auto radius =
                            math.add(math.abs(math.sub(output, centre(approximated, false))),
                                     centre(approximated, true));
                        if (math.compare(radius, maximum_error) > 0)
                            maximum_error = std::move(radius);
                    }
                }
            }
            auto checking = limits.chart_verification;
            checking.cancellation = limits.arithmetic.cancellation;
            math.reserve_nested_rationals(checking.algebra, partition_add(growing, 1));
            auto result = partition_take(chart.lift_rhs(
                {identity, independent}, {identity, rhs.values()}, checking, workspace, plan));
            return CaseStructuralExactLift{std::move(rhs), std::move(result),
                                           std::move(maximum_error)};
        };
        try {
            return compute();
        } catch (const m::ArithmeticFailure& failure) {
            return std::unexpected(core::Failure{failure.code, operation});
        }
    });
}

} // namespace uvlas::coupling::detail
