#include "case_pg_chart.hpp"

#include "case_pg_support_input.hpp"
#include "structure/pg_reference_domain.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <memory>
#include <memory_resource>
#include <type_traits>

namespace uvlas::coupling::detail {
namespace {
namespace pg = structure;
namespace s = io::schema;
using core::FailureCode;
using namespace pg::detail;
using R = pg::ExactRational;
constexpr auto operation = core::Operation::PrepareModel;

struct Input {
    const io::SiCase& source;
    ReferenceArithmetic math;
    std::uint64_t work;

    void step(std::uint64_t count = 1) {
        math.cancel();
        if (count > work)
            partition_fail(FailureCode::ResourceExceeded);
        work -= count;
    }

    template <class T>
    auto allocate(std::uint64_t count) {
        step(count);
        return partition_take(math.plan.allocate<T>(partition_index(count), operation));
    }

    template <class Field>
    R rational(Field input) {
        const auto value = input.record();
        const auto numerator = value.numerator.text();
        const auto denominator = value.denominator.text();
        step(partition_add(1, partition_add(numerator.size(), denominator.size())));
        math.charge();
        return partition_take(
            R::from_integer_pair(numerator, denominator, math.leaf_limits, math.plan));
    }

    template <class T>
    auto views(std::uint64_t count) {
        static_assert(std::is_trivially_destructible_v<T>);
        return allocate<std::byte>(partition_add(partition_mul(count, sizeof(T)), alignof(T) - 1));
    }
};

pg::PgPreparedPatch prepare_patch(const io::SiCase& source, const CasePgRuleBinding& rule,
                                  const CasePgChartLimits& limits, core::Workspace& workspace,
                                  core::ResourcePlan& plan) {
    auto geometry = partition_take(
        prepare_case_pg_reference(source, rule.shell, limits.reference, workspace, plan));
    const auto owners = partition_take(geometry.source());
    auto rank = partition_take(pg::PgRankPreparation::prepare(owners.dictionary, owners.partition,
                                                              limits.rank, workspace, plan));
    auto basis = partition_take(pg::PgReferenceBasis::prepare(std::move(geometry), std::move(rank),
                                                              limits.rank, workspace, plan));
    auto gram = partition_take(pg::PgDiscreteGram::prepare(
        std::move(basis), rule.points, limits.gram_grid, limits.gram, workspace, plan));
    auto certificate = partition_take(
        std::move(gram).generate_and_certify(limits.factor_grid, limits.factor, workspace, plan));
    return partition_take(
        pg::PgPreparedPatch::prepare(std::move(certificate), limits.patch, workspace, plan));
}

struct Endpoint {
    core::Index patch{};
    core::Index boundary{};
};

Endpoint endpoint(const s::Region& region, const s::Model& model, Input& input) {
    const auto target = region.target.record();
    if (target.kind.selected() != decltype(target.kind.selected())::ItemShell ||
        region.kind.selected() != decltype(region.kind.selected())::ItemEdge ||
        region.side.selected() != decltype(region.side.selected())::ItemInterior)
        partition_fail(FailureCode::Unsupported);
    for (core::Index i = 0; i < partition_index(model.shells.size()); ++i) {
        input.step();
        const auto shell = partition_take(model.shells.at(i)).record();
        if (shell.id.unsigned_value() != target.id.unsigned_value())
            continue;
        if (shell.geometry.unsigned_value() != region.geometry.unsigned_value())
            partition_fail(FailureCode::InvalidInput);
        return {i, partition_index(region.boundary.unsigned_value())};
    }
    partition_fail(FailureCode::InvalidInput);
}

struct NamedBoundary {
    core::Index ordinal;
    std::span<const pg::PgPartitionPoint> vertices;
    std::span<const pg::PgPartitionEdge> edges;
    std::span<const pg::PgPartitionPiece> pieces;
    std::span<const pg::PgPartitionEdgeTrace> traces;
    std::span<const pg::PgPartitionSource> sources;
    std::span<const pg::PgPartitionPoint> original_points;
    pg::PgPartitionBoundaryRecord record;
};

NamedBoundary named_boundary(const pg::PgPreparedPatch& patch, core::Index ordinal) {
    const auto basis = partition_take(patch.conditioning().reference_basis());
    const auto& partition = partition_take(basis->source()).partition;
    return {ordinal,
            partition_take(partition.vertices()),
            partition_take(partition.edges()),
            partition_take(partition.pieces()),
            partition_take(partition.edge_traces()),
            partition_take(partition.sources()),
            partition_take(partition.source_points()),
            reference_at(partition_take(partition.boundaries()), ordinal)};
}

const pg::PgPartitionNamedSource* named_segment(const NamedBoundary& boundary,
                                                const pg::PgPartitionEdgeTrace& trace) {
    const auto& descriptor = reference_at(boundary.sources, trace.source).descriptor;
    const auto* segment = std::get_if<pg::PgPartitionNamedSource>(&descriptor);
    return segment != nullptr && segment->boundary == boundary.ordinal ? segment : nullptr;
}

R index_rational(core::Index value, ReferenceArithmetic& math) {
    std::array<char, 24> digits{};
    const auto result = std::to_chars(digits.data(), digits.data() + digits.size(), value);
    return math.constant({digits.data(), static_cast<std::size_t>(result.ptr - digits.data())});
}

std::uint64_t named_size(const NamedBoundary& boundary, Input& input) {
    std::uint64_t count = 0;
    for (const auto& trace : boundary.traces) {
        input.step();
        if (named_segment(boundary, trace) != nullptr)
            ++count;
    }
    if (count == 0)
        partition_fail(FailureCode::GeometryInvalid);
    return count;
}

void admit_trace_input(Input& input, std::uint64_t live) {
    auto& math = input.math;
    // The preceding constraint's staging buffers are already gone.
    math.leaf_limits = math.limits.algebra;
    math.admit(partition_add(reference_live_rationals, live));
}

struct EdgeRange {
    core::Index edge{};
    pg::PgCellIndex cell;
    R first, last;
    bool increasing{}, forward{};
};

template <class Parameter>
auto ranges(const NamedBoundary& boundary, const std::uint64_t count, Parameter parameter,
            Input& input) {
    auto records = input.allocate<EdgeRange>(count);
    std::size_t n = 0;
    auto& math = input.math;
    for (const auto& trace : boundary.traces) {
        input.step();
        const auto* segment = named_segment(boundary, trace);
        if (segment == nullptr)
            continue;
        const auto& edge = reference_at(boundary.edges, trace.edge);
        const auto& a = reference_at(boundary.vertices, edge.vertices[0]);
        const auto& b = reference_at(boundary.vertices, edge.vertices[1]);
        R x = parameter(*segment, a);
        R y = parameter(*segment, b);
        const auto order = math.compare(x, y);
        if (order == 0)
            partition_fail(FailureCode::GeometryInvalid);
        const bool increasing = order < 0;
        if (!increasing)
            std::swap(x, y);
        records.values()[n++] = {trace.edge,   reference_at(boundary.pieces, trace.piece).cell,
                                 std::move(x), std::move(y),
                                 increasing,   increasing == trace.canonical_direction};
    }
    std::ranges::sort(records.values(), [&](const auto& a, const auto& b) {
        input.step();
        return math.compare(a.first, b.first) < 0;
    });
    for (std::size_t i = 1; i < records.values().size(); ++i)
        if (math.compare(records.values()[i - 1].last, records.values()[i].first) != 0)
            partition_fail(FailureCode::GeometryInvalid);
    return records;
}

std::span<const pg::PgTraceEdgeSpan> construct_spans(std::span<const EdgeRange> records,
                                                     std::span<std::byte> storage, const R& zero,
                                                     const R& one) {
    std::pmr::monotonic_buffer_resource arena(storage.data(), storage.size(),
                                              std::pmr::null_memory_resource());
    auto* result =
        std::pmr::polymorphic_allocator<pg::PgTraceEdgeSpan>{&arena}.allocate(records.size());
    for (std::size_t i = 0; i < records.size(); ++i) {
        const auto& record = records[i];
        std::construct_at(result + i, record.edge, record.first, record.last,
                          record.increasing ? zero : one, record.increasing ? one : zero,
                          std::nullopt);
    }
    return {result, records.size()};
}

struct Affine {
    std::array<R, 2> value, slope;
};

template <class Map>
Affine read_map(Map map, const R& one, Input& input) {
    if (map.size() != 3)
        partition_fail(FailureCode::InvalidInput);
    Affine output;
    for (core::Index d = 0; d < 3; ++d) {
        const auto terms = partition_take(map.at(d));
        if (terms.size() > 2 || (d == 2 && terms.size() != 1))
            partition_fail(FailureCode::Unsupported);
        std::array<bool, 2> seen{};
        for (core::Index i = 0; i < partition_index(terms.size()); ++i) {
            input.step();
            const auto term = partition_take(terms.at(i)).record();
            const auto eta = partition_take(term.powers.at(0)).unsigned_value();
            const auto zeta = partition_take(term.powers.at(1)).unsigned_value();
            auto value = input.rational(term.coefficient);
            if (d == 2) {
                if (eta != 0 || zeta != 1 || input.math.compare(value, one) != 0)
                    partition_fail(FailureCode::Unsupported);
                continue;
            }
            if (eta > 1 || zeta != 0 || seen[eta])
                partition_fail(FailureCode::Unsupported);
            seen[eta] = true;
            const auto k = static_cast<std::size_t>(d);
            (eta == 0 ? output.value[k] : output.slope[k]) = std::move(value);
        }
    }
    if (output.slope[0].sign() == 0 && output.slope[1].sign() == 0)
        partition_fail(FailureCode::GeometryInvalid);
    return output;
}

bool same_map(const Affine& a, const Affine& b, ReferenceArithmetic& math) {
    for (std::size_t d = 0; d < a.value.size(); ++d)
        if (math.compare(a.value[d], b.value[d]) != 0 || math.compare(a.slope[d], b.slope[d]) != 0)
            return false;
    return true;
}

R parameter(const Affine& line, const pg::PgPartitionPoint& point, ReferenceArithmetic& math) {
    const auto d = line.slope[0].sign() == 0 ? 1U : 0U;
    const auto& first = d == 0 ? point.s : point.t;
    const auto& second = d == 0 ? point.t : point.s;
    auto eta = math.div(math.sub(first, line.value[d]), line.slope[d]);
    const auto expected = math.add(line.value[1 - d], math.mul(eta, line.slope[1 - d]));
    if (math.compare(second, expected) != 0)
        partition_fail(FailureCode::GeometryInvalid);
    return eta;
}

struct Tile {
    R first, last;
    pg::PgCellIndex left{}, right{};
    bool reversed{};
};

Tile tile(const s::TraceMap& mapping, const R& one, Input& input) {
    auto& math = input.math;
    if (mapping.domain.size() != 4)
        partition_fail(FailureCode::Unsupported);
    Tile result;
    const auto negative = math.neg(one);
    std::array<R, 4> eta;
    std::array<bool, 4> upper{};
    for (core::Index i = 0; i < 4; ++i) {
        const auto point = partition_take(mapping.domain.at(i));
        const auto j = static_cast<std::size_t>(i);
        eta[j] = input.rational(partition_take(point.at(0)));
        const auto depth = input.rational(partition_take(point.at(1)));
        upper[j] = math.compare(depth, one) == 0;
        if (!upper[j] && math.compare(depth, negative) != 0)
            partition_fail(FailureCode::Unsupported);
    }
    const auto bounds = std::ranges::minmax_element(
        eta, [&](const auto& a, const auto& b) { return math.compare(a, b) < 0; });
    result.first = math.copy(*bounds.min);
    result.last = math.copy(*bounds.max);
    if (math.compare(result.first, result.last) == 0)
        partition_fail(FailureCode::GeometryInvalid);
    std::array<bool, 4> seen{};
    for (std::size_t i = 0; i < eta.size(); ++i) {
        const bool right = math.compare(eta[i], result.last) == 0;
        if (!right && math.compare(eta[i], result.first) != 0)
            partition_fail(FailureCode::Unsupported);
        const auto corner = (right ? 2U : 0U) + (upper[i] ? 1U : 0U);
        if (seen[corner])
            partition_fail(FailureCode::InvalidInput);
        seen[corner] = true;
        auto next = (i + 1) % eta.size();
        if ((upper[i] == upper[next]) == (math.compare(eta[i], eta[next]) == 0))
            partition_fail(FailureCode::GeometryInvalid);
    }
    const auto cell = [](auto field) {
        return pg::PgCellIndex{partition_index(partition_take(field.at(0)).unsigned_value()),
                               partition_index(partition_take(field.at(1)).unsigned_value())};
    };
    result.left = cell(mapping.first_cell);
    result.right = cell(mapping.second_cell);
    result.reversed =
        mapping.orientation.selected() == decltype(mapping.orientation.selected())::ItemReversed;
    return result;
}

void validate_coverage(std::span<const EdgeRange> ranges_, std::span<const Tile> tiles,
                       bool first_side, bool relative_reversed, Input& input) {
    auto& math = input.math;
    if (math.compare(ranges_.front().first, tiles.front().first) != 0 ||
        math.compare(ranges_.back().last, tiles.back().last) != 0)
        partition_fail(FailureCode::GeometryInvalid);
    for (std::size_t i = 0, j = 0; i < ranges_.size() && j < tiles.size();) {
        input.step();
        const auto& range = ranges_[i];
        const auto& tile_ = tiles[j];
        if (range.cell != (first_side ? tile_.left : tile_.right) ||
            range.forward != ranges_.front().forward || tile_.reversed != relative_reversed)
            partition_fail(FailureCode::GeometryInvalid);
        const auto end = math.compare(range.last, tile_.last);
        if (end <= 0)
            ++i;
        if (end >= 0)
            ++j;
    }
}

pg::PgFullTrace connection(const s::JointAdhesive& joint, const pg::PgPreparedPatch& left,
                           core::Index left_boundary, const pg::PgPreparedPatch& right,
                           core::Index right_boundary, const pg::PgTraceLimits& trace_limits,
                           core::Workspace& workspace, Input& input) {
    if (joint.split_at_attachment.boolean() || joint.mapping.size() == 0)
        partition_fail(FailureCode::Unsupported);
    auto& math = input.math;
    const auto left_named = named_boundary(left, left_boundary);
    const auto right_named = named_boundary(right, right_boundary);
    const auto lcount = named_size(left_named, input);
    const auto rcount = named_size(right_named, input);
    const auto records = partition_add(partition_add(lcount, rcount), joint.mapping.size());
    admit_trace_input(input, partition_mul(2, records));
    const auto zero = R{};
    const auto one = math.constant("1");
    const auto first = partition_take(joint.mapping.at(0)).record();
    const auto left_map = read_map(first.first, one, input);
    const auto right_map = read_map(first.second, one, input);
    auto tiles = input.allocate<Tile>(joint.mapping.size());
    for (std::size_t i = 0; i < tiles.values().size(); ++i) {
        input.step();
        const auto mapping = partition_take(joint.mapping.at(partition_index(i))).record();
        if (!same_map(left_map, read_map(mapping.first, one, input), math) ||
            !same_map(right_map, read_map(mapping.second, one, input), math))
            partition_fail(FailureCode::Unsupported);
        tiles.values()[i] = tile(mapping, one, input);
    }
    std::ranges::sort(tiles.values(), [&](const auto& a, const auto& b) {
        input.step();
        return math.compare(a.first, b.first) < 0;
    });
    for (std::size_t i = 1; i < tiles.values().size(); ++i)
        if (math.compare(tiles.values()[i - 1].last, tiles.values()[i].first) != 0)
            partition_fail(FailureCode::GeometryInvalid);

    const auto l = ranges(
        left_named, lcount,
        [&](const auto&, const auto& point) { return parameter(left_map, point, math); }, input);
    const auto r = ranges(
        right_named, rcount,
        [&](const auto&, const auto& point) { return parameter(right_map, point, math); }, input);
    const auto reversed = l.values().front().forward != r.values().front().forward;
    validate_coverage(l.values(), tiles.values(), true, reversed, input);
    validate_coverage(r.values(), tiles.values(), false, reversed, input);
    auto lviews = input.views<pg::PgTraceEdgeSpan>(l.values().size());
    auto rviews = input.views<pg::PgTraceEdgeSpan>(r.values().size());
    const auto lspans = construct_spans(l.values(), lviews.values(), zero, one);
    const auto rspans = construct_spans(r.values(), rviews.values(), zero, one);
    auto nested = trace_limits;
    math.reserve_nested_rationals(nested.arithmetic.algebra,
                                  partition_add(10, partition_mul(2, records)));
    auto trace = partition_take(
        pg::PgFullTrace::prepare(left, lspans, right, rspans, nested, workspace, math.plan));
    const auto& overlaps = trace.overlaps();
    // The declared material map pairs equal zeta. It can represent equal depths
    // only when both physical full thicknesses match across the entire edge.
    for (const auto& overlap : overlaps) {
        input.step();
        const auto& a = overlap.left_depth;
        const auto& b = overlap.right_depth;
        if (math.compare(a.lower_first, b.lower_first) != 0 ||
            math.compare(a.lower_last, b.lower_last) != 0 ||
            math.compare(a.upper_first, b.upper_first) != 0 ||
            math.compare(a.upper_last, b.upper_last) != 0)
            partition_fail(FailureCode::GeometryInvalid);
    }
    return trace;
}

bool validate_support(const s::SupportTraceConstraint& declaration,
                      const CaseMotionFrameLimits& limits, bool motion_enabled, Input& input) {
    if (declaration.components.size() != 3)
        partition_fail(FailureCode::Unsupported);
    input.step(7); // The event count and at most six component lookups.
    const auto stationary =
        partition_take(stationary_case_pg_motion(input.source, declaration.motion.record()));
    if (!motion_enabled && !stationary)
        partition_fail(FailureCode::Unsupported);
    std::array<bool, 3> seen{};
    for (core::Index i = 0; i < 3; ++i) {
        input.step();
        const auto d = partition_take(declaration.components.at(i)).unsigned_value();
        if (d >= seen.size() || seen[d])
            partition_fail(FailureCode::InvalidInput);
        seen[d] = true;
    }
    if (!motion_enabled) {
        const auto frame = partition_take(CaseMotionFrame::prepare(input.source, declaration.motion,
                                                                   limits, input.math.plan,
                                                                   input.math.limits.cancellation));
        static_cast<void>(frame);
    }
    return stationary;
}

pg::PgFullTrace support(const pg::PgPreparedPatch& patch, core::Index ordinal,
                        const pg::PgTraceLimits& trace_limits, core::Workspace& workspace,
                        Input& input) {
    const auto boundary = named_boundary(patch, ordinal);
    const auto end = boundary.record.point_count - 1;
    auto& math = input.math;
    const auto count = named_size(boundary, input);
    admit_trace_input(input, partition_mul(2, count));
    const auto records = ranges(
        boundary, count,
        [&](const pg::PgPartitionNamedSource& segment, const pg::PgPartitionPoint& point) {
            const auto i = boundary.record.first_point + segment.segment;
            const auto& a = reference_at(boundary.original_points, i);
            const auto& b = reference_at(boundary.original_points, i + 1);
            const bool along_s = math.compare(a.s, b.s) != 0;
            const auto& first = along_s ? a.s : a.t;
            const auto& last = along_s ? b.s : b.t;
            const auto& value = along_s ? point.s : point.t;
            const auto local = math.div(math.sub(value, first), math.sub(last, first));
            return math.add(index_rational(segment.segment, math), local);
        },
        input);
    const auto zero = R{};
    const auto one = math.constant("1");
    if (math.compare(records.values().front().first, zero) != 0 ||
        math.compare(records.values().back().last, index_rational(end, math)) != 0)
        partition_fail(FailureCode::GeometryInvalid);
    auto views = input.views<pg::PgTraceEdgeSpan>(records.values().size());
    const auto spans = construct_spans(records.values(), views.values(), zero, one);
    auto nested = trace_limits;
    math.reserve_nested_rationals(nested.arithmetic.algebra,
                                  partition_add(2, partition_mul(2, count)));
    return partition_take(pg::PgFullTrace::prepare_support(
        patch, spans, partition_take(patch.info()).identity, patch.reference_coordinates(), nested,
        workspace, math.plan));
}

struct MotionPreparation {
    aero::detail::pchip::PchipInputStamp stamp;
    const CaseStructuralSupportLimits& limits;
    aero::detail::pchip::CurveWorkspace& workspace;
};

} // namespace

core::Result<const pg::PgPreparedPatch*> CasePgChart::patch(core::Index index) const& noexcept {
    if (index < 0 || static_cast<std::uint64_t>(index) >= patches_.values().size())
        return std::unexpected(core::Failure{FailureCode::InvalidInput, operation});
    const auto& patch = patches_.values()[static_cast<std::size_t>(index)];
    if (!patch)
        return std::unexpected(core::Failure{FailureCode::InvalidInput, operation});
    return &*patch;
}

const CaseStructuralSupport* CasePgChart::motion() const& noexcept {
    return std::get_if<CaseStructuralSupport>(&chart_);
}

const pg::PgSparseChart& CasePgChart::chart() const& noexcept {
    if (const auto* ready = motion())
        return ready->support().chart();
    return std::get<pg::PgSparseChart>(chart_);
}

namespace case_pg_chart_detail {
struct Builder {
    static core::Result<CasePgChart> prepare(const io::SiCase&, std::span<const CasePgRuleBinding>,
                                             const CasePgChartLimits&, const MotionPreparation*,
                                             core::Workspace&, core::ResourcePlan&,
                                             std::stop_token) noexcept;
};

core::Result<CasePgChart> Builder::prepare(const io::SiCase& source,
                                           std::span<const CasePgRuleBinding> rules,
                                           const CasePgChartLimits& caller,
                                           const MotionPreparation* motion,
                                           core::Workspace& workspace, core::ResourcePlan& plan,
                                           std::stop_token stop) noexcept {
    return core::guard_entry(
        operation, FailureCode::InvalidInput, {}, [&]() -> core::Result<CasePgChart> {
            try {
                std::stop_source cancellation;
                const auto request = [&]() noexcept { cancellation.request_stop(); };
                const std::stop_callback caller_stop{stop, request};
                const std::stop_callback reference_stop{caller.reference.cancellation, request};
                const std::stop_callback graph_stop{caller.reference.graph.cancellation, request};
                const std::stop_callback coefficient_stop{
                    caller.reference.coefficients.cancellation, request};
                const std::stop_callback partition_stop{caller.reference.partition.cancellation,
                                                        request};
                const std::stop_callback geometry_stop{caller.reference.geometry.cancellation,
                                                       request};
                const std::stop_callback rank_stop{caller.rank.cancellation, request};
                const std::stop_callback gram_stop{caller.gram.arithmetic.cancellation, request};
                const std::stop_callback factor_stop{caller.factor.cancellation, request};
                const std::stop_callback patch_stop{caller.patch.cancellation, request};
                const std::stop_callback trace_stop{caller.trace.arithmetic.cancellation, request};
                const std::stop_callback constraints_stop{
                    caller.constraints.arithmetic.cancellation, request};
                const std::stop_callback chart_stop{caller.chart.arithmetic.cancellation, request};

                const std::stop_callback motion_trace_stop{
                    motion == nullptr ? std::stop_token{}
                                      : motion->limits.region_trace.arithmetic.cancellation,
                    request};
                const std::stop_callback motion_exact_stop{
                    motion == nullptr ? std::stop_token{}
                                      : motion->limits.exact.arithmetic.cancellation,
                    request};
                const std::stop_callback motion_chart_stop{
                    motion == nullptr ? std::stop_token{}
                                      : motion->limits.exact.chart_verification.cancellation,
                    request};

                const auto token = cancellation.get_token();
                auto limits = caller;
                limits.reference.cancellation = token;
                limits.rank.cancellation = token;
                limits.gram.arithmetic.cancellation = token;
                limits.factor.cancellation = token;
                limits.patch.cancellation = token;
                limits.trace.arithmetic.cancellation = token;
                limits.constraints.arithmetic.cancellation = token;
                limits.chart.arithmetic.cancellation = token;

                Input input{source, {limits.trace.arithmetic, plan}, limits.maximum_input_work};
                input.step();
                const auto model = partition_take(partition_take(source.source_case()).get().root())
                                       .model.record();
                const auto source_count = partition_add(model.joints.size(), model.supports.size());
                if (rules.empty() || rules.size() != model.shells.size())
                    partition_fail(FailureCode::InvalidInput);
                if (rules.size() > limits.constraints.maximum_patches ||
                    source_count > limits.constraints.maximum_sources)
                    partition_fail(FailureCode::ResourceExceeded);

                for (std::size_t i = 0; i < rules.size(); ++i) {
                    input.step();
                    const auto shell = partition_take(model.shells.at(partition_index(i))).record();
                    if (rules[i].shell != shell.id.unsigned_value() || rules[i].points.empty())
                        partition_fail(FailureCode::InvalidInput);
                }
                const auto support_count = static_cast<std::size_t>(model.supports.size());
                core::PlannedBuffer<std::byte> support_storage;
                if (motion != nullptr) {
                    if (support_count > motion->limits.maximum_group_supports)
                        partition_fail(FailureCode::ResourceExceeded);
                    support_storage = input.views<pg::PgConstraintId>(support_count);
                }
                std::pmr::monotonic_buffer_resource support_arena(support_storage.values().data(),
                                                                  support_storage.values().size(),
                                                                  std::pmr::null_memory_resource());
                auto* support_ids =
                    motion == nullptr
                        ? nullptr
                        : std::pmr::polymorphic_allocator<pg::PgConstraintId>{&support_arena}
                              .allocate(support_count);
                std::size_t moving_count = 0;
                std::size_t stationary_start = support_count;
                for (core::Index i = 0; i < partition_index(model.supports.size()); ++i) {
                    input.step();
                    const auto constraint = partition_take(model.supports.at(i)).record();
                    const auto* fixed = std::get_if<s::SupportTraceConstraint>(&constraint);
                    if (fixed == nullptr)
                        partition_fail(FailureCode::Unsupported);
                    static_cast<void>(endpoint(fixed->region.record(), model, input));
                    const bool stationary =
                        validate_support(*fixed, limits.support_frame, motion != nullptr, input);
                    if (support_ids != nullptr) {
                        const auto id = pg::PgConstraintId{fixed->id.unsigned_value()};
                        const auto offset = stationary ? --stationary_start : moving_count++;
                        std::construct_at(support_ids + offset, id);
                    }
                }
                if (motion != nullptr && moving_count == 0)
                    partition_fail(FailureCode::InvalidInput);
                auto patches = input.allocate<std::optional<pg::PgPreparedPatch>>(rules.size());
                auto traces = input.allocate<std::optional<pg::PgFullTrace>>(source_count);
                auto patch_storage = input.views<pg::PgConstraintPatchBinding>(rules.size());
                auto trace_storage = input.views<pg::PgConstraintTraceBinding>(source_count);
                std::pmr::monotonic_buffer_resource patch_arena(patch_storage.values().data(),
                                                                patch_storage.values().size(),
                                                                std::pmr::null_memory_resource());
                std::pmr::monotonic_buffer_resource trace_arena(trace_storage.values().data(),
                                                                trace_storage.values().size(),
                                                                std::pmr::null_memory_resource());
                auto* patch_bindings =
                    std::pmr::polymorphic_allocator<pg::PgConstraintPatchBinding>{&patch_arena}
                        .allocate(rules.size());
                auto* trace_bindings =
                    std::pmr::polymorphic_allocator<pg::PgConstraintTraceBinding>{&trace_arena}
                        .allocate(static_cast<std::size_t>(source_count));

                for (std::size_t i = 0; i < rules.size(); ++i) {
                    input.step();
                    const auto& patch = patches.values()[i].emplace(
                        prepare_patch(source, rules[i], limits, workspace, plan));
                    std::construct_at(patch_bindings + i, numerics::PatchId{rules[i].shell}, patch);
                }

                std::size_t k = 0;
                for (core::Index i = 0; i < partition_index(model.joints.size()); ++i, ++k) {
                    input.step();
                    const auto joint = partition_take(model.joints.at(i)).record();
                    const auto* adhesive = std::get_if<s::JointAdhesive>(&joint);
                    if (adhesive == nullptr)
                        partition_fail(FailureCode::Unsupported);
                    const auto l = endpoint(adhesive->first.record(), model, input);
                    const auto r = endpoint(adhesive->second.record(), model, input);
                    if (l.patch == r.patch)
                        partition_fail(FailureCode::InvalidInput);
                    const auto& trace = traces.values()[k].emplace(connection(
                        *adhesive, patch_bindings[static_cast<std::size_t>(l.patch)].patch,
                        l.boundary, patch_bindings[static_cast<std::size_t>(r.patch)].patch,
                        r.boundary, limits.trace, workspace, input));
                    std::construct_at(trace_bindings + k,
                                      pg::PgConstraintId{adhesive->id.unsigned_value()}, trace,
                                      l.patch, r.patch);
                }
                for (core::Index i = 0; i < partition_index(model.supports.size()); ++i, ++k) {
                    input.step();
                    const auto constraint = partition_take(model.supports.at(i)).record();
                    const auto* fixed = std::get_if<s::SupportTraceConstraint>(&constraint);
                    if (fixed == nullptr)
                        partition_fail(FailureCode::Unsupported);
                    const auto selected = endpoint(fixed->region.record(), model, input);
                    const auto& trace = traces.values()[k].emplace(
                        support(patch_bindings[static_cast<std::size_t>(selected.patch)].patch,
                                selected.boundary, limits.trace, workspace, input));
                    std::construct_at(trace_bindings + k,
                                      pg::PgConstraintId{fixed->id.unsigned_value()}, trace,
                                      selected.patch, std::nullopt);
                }
                auto original = partition_take(pg::PgGlobalConstraints::prepare(
                    {patch_bindings, rules.size()}, {trace_bindings, k}, limits.constraints,
                    workspace, plan));
                auto chart = partition_take(
                    pg::PgSparseChart::prepare(std::move(original), limits.chart, workspace, plan));
                auto result = CasePgChart{std::move(patches), std::move(chart)};
                if (motion != nullptr) {
                    // Bind the very chart and patch owners used in row assembly. No new trace
                    // or factorization replaces this chart during evaluation.
                    const CasePgSupportSelection selection{
                        {support_ids, moving_count},
                        {support_ids + stationary_start, support_count - stationary_start}};
                    auto support = partition_take(prepare_case_structural_pg_support(
                        source, std::move(std::get<pg::PgSparseChart>(result.chart_)),
                        {patch_bindings, rules.size()}, selection, motion->stamp, workspace,
                        motion->workspace, motion->limits, plan, token));
                    result.chart_.emplace<CaseStructuralSupport>(std::move(support));
                }
                input.step(0);
                return result;
            } catch (const pg::detail::PartitionAbort& error) {
                return std::unexpected(error.failure);
            }
        });
}

} // namespace case_pg_chart_detail

core::Result<CasePgChart>
prepare_case_pg_chart(const io::SiCase& source, std::span<const CasePgRuleBinding> rules,
                      const CasePgChartLimits& limits, core::Workspace& workspace,
                      core::ResourcePlan& plan, std::stop_token stop) noexcept {
    return case_pg_chart_detail::Builder::prepare(source, rules, limits, nullptr, workspace, plan,
                                                  stop);
}

core::Result<CasePgChart> prepare_case_pg_motion_chart(
    const io::SiCase& source, std::span<const CasePgRuleBinding> rules,
    const CasePgChartLimits& limits, aero::detail::pchip::PchipInputStamp stamp,
    const CaseStructuralSupportLimits& support_limits, core::Workspace& workspace,
    aero::detail::pchip::CurveWorkspace& motion_workspace, core::ResourcePlan& plan,
    std::stop_token stop) noexcept {
    const MotionPreparation motion{stamp, support_limits, motion_workspace};
    return case_pg_chart_detail::Builder::prepare(source, rules, limits, &motion, workspace, plan,
                                                  stop);
}

} // namespace uvlas::coupling::detail
