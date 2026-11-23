#include "p1_induction.hpp"

#include "geometry_snapshot_storage.hpp"
#include "local_geometry_bounds_private.hpp"
#include "p1_exact_product_private.hpp"
#include "p1_gauss_private.hpp"
#include "p1_program_private.hpp"
#include "p1_quadrature_private.hpp"
#include "p1_qualification_private.hpp"
#include "p1_shared_kernel.hpp"

#include <bit>
#include <cfenv>
#include <cmath>
#include <new>

namespace uvlas::aero::detail {
namespace {
constexpr auto operation = core::Operation::EvaluateAero;
bool supported_environment() noexcept {
    if (std::fegetround() != FE_TONEAREST)
        return false;
    // Both producing and consuming subnormals must retain their exact bits.
    const volatile double normal = std::numeric_limits<double>::min();
    const volatile double half = 0.5;
    const volatile double subnormal = std::numeric_limits<double>::denorm_min();
    const volatile double one = 1.0;
    const double produced = normal * half;
    const double consumed = subnormal * one;
    return std::bit_cast<std::uint64_t>(produced) == UINT64_C(0x0008000000000000) &&
           std::bit_cast<std::uint64_t>(consumed) == 1;
}
constexpr std::uint64_t input_work = 13 * 4 + 4 * 3 * 4 + 3 * 3 * 6 + 3 * 3 * 6 + 3 * 3 + 3 * 3;
constexpr std::uint64_t table_work = 8 * (5 * 8 + 8) + 12 * (5 * 12 + 8);
constexpr std::uint64_t constant_work = 781 + 393;
constexpr std::uint64_t result_work = UINT64_C(126) * 14;
struct Shape {
    std::size_t nodes = 0;
    std::uint64_t work = 0;
};
Shape maximum_shape() {
    Shape maximum;
    for (unsigned pattern = 0; pattern < 729U; ++pattern) {
        unsigned remaining = pattern;
        p1::Options options;
        for (std::size_t i = 0; i < 3; ++i) {
            options.start_sign[i] = static_cast<std::int8_t>(static_cast<int>(remaining % 3U) - 1);
            remaining /= 3U;
            options.end_sign[i] = static_cast<std::int8_t>(static_cast<int>(remaining % 3U) - 1);
            remaining /= 3U;
        }
        for (const bool exterior : {false, true}) {
            options.exterior = exterior;
            p1::Program count{{}};
            const auto output = p1::near_program(count, options);
            (void)output;
            maximum.nodes = std::max(maximum.nodes, count.count);
            maximum.work = std::max(maximum.work, count.work);
        }
    }
    p1::Program far{{}};
    const auto output = p1::far_program(far);
    (void)output;
    const auto far_nodes = p1::checked_sum(far.count, 18);
    if (far_nodes > static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max()))
        throw p1::ArithmeticFailure{core::FailureCode::ResourceExceeded};
    maximum.nodes = std::max(maximum.nodes, static_cast<std::size_t>(far_nodes));
    constexpr std::uint64_t scalar_root = p1::root_steps;
    // Every sample weights once and adds all 9x14 returned jet components.
    constexpr std::uint64_t accumulation = UINT64_C(208) * (1 + UINT64_C(9) * 14);
    // Each normalized vector uses 3x7 scalar operations, two norm roots,
    // and one relative-budget multiplication. ell and ell^2 are prepared once.
    constexpr std::uint64_t comparison =
        42 * (UINT64_C(3) * 7 + 2 * scalar_root + 1) + scalar_root + 1;
    maximum.work =
        std::max(maximum.work, p1::checked_sum(p1::checked_product(208, far.work),
                                               p1::checked_sum(accumulation, comparison)));
    maximum.work =
        p1::checked_sum(maximum.work, input_work + table_work + constant_work + result_work);
    return maximum;
}
OwnedGeometryFailure failed(core::FailureCode code, const GeometryContext& context) noexcept {
    auto details = context.diagnostics;
    if (!details.model)
        details.model = context.model;
    if (!details.numerical_generation)
        details.numerical_generation = context.numerical;
    return {core::Failure{code, operation, details}, {}, context};
}
core::Failure with_details(core::Failure original, const core::FailureDetails& supplied) noexcept {
    auto& target = original.details;
    if (!target.model)
        target.model = supplied.model;
    if (!target.parent)
        target.parent = supplied.parent;
    if (!target.trial)
        target.trial = supplied.trial;
    if (!target.numerical_generation)
        target.numerical_generation = supplied.numerical_generation;
    if (!target.component)
        target.component = supplied.component;
    if (!target.problem)
        target.problem = supplied.problem;
    if (!target.resources)
        target.resources = supplied.resources;
    if (!target.last_accepted)
        target.last_accepted = supplied.last_accepted;
    return original;
}
template <class T>
T needed(AeroResult<T> result) {
    if (!result)
        throw result.error();
    return std::move(*result);
}
void decode_words(const ScalarWords& words, ExactInteger& out, ExactInteger& temporary,
                  p1::Work& work) {
    work.step();
    if (const auto* one = std::get_if<SingleWord>(&words))
        local_decode(one->value, out);
    else {
        const auto& pair = std::get<DoubleWord>(words);
        local_decode(pair.high, out);
        local_decode(pair.low, temporary);
        out += temporary;
    }
}
std::int8_t sign(const ExactInteger& value) {
    return value == 0 ? std::int8_t{0} : (value < 0 ? std::int8_t{-1} : std::int8_t{1});
}
p1::Interval word(double value) {
    ExactInteger decoded;
    local_decode(value, decoded);
    return p1::exact(decoded, -1074);
}
using Rule = p1::QuadratureRule;
} // namespace
struct P1Scratch {
    Shape shape;
    core::PlannedBuffer<p1::Node> program;
    core::PlannedBuffer<p1::Jet> jets;
    std::array<std::array<ExactInteger, 3>, 4> points;
    std::array<ExactInteger, 3> weights;
    std::array<ExactInteger, 6> temporary;
    LocalRootScratch output;
    p1::ExactProductScratch product;
    Rule low, high;
    std::array<p1::Interval, 16> input;
    p1::FarComponents low_values, high_values;
    p1::Interval pi, ln2;
};
namespace {
std::uint64_t shape_payload(const Shape& shape) {
    if (shape.nodes > static_cast<std::size_t>(std::numeric_limits<core::Index>::max()))
        throw p1::ArithmeticFailure{core::FailureCode::ResourceExceeded};
    const auto count = static_cast<core::Index>(shape.nodes);
    const auto nodes = core::checked_byte_count(count, sizeof(p1::Node), operation);
    const auto jets = core::checked_byte_count(count, sizeof(p1::Jet), operation);
    if (!nodes || !jets)
        throw p1::ArithmeticFailure{core::FailureCode::ResourceExceeded};
    return p1::checked_sum(sizeof(P1Scratch), p1::checked_sum(*nodes, *jets));
}
} // namespace
struct P1Data {
    OwnedGeometry geometry;
    QueryRequest query;
    LocalAmbientDescriptor ambient;
    P1Accuracy accuracy;
    P1Branch branch;
    P1EvaluationDiagnostic diagnostic;
    P1WorkPlan plan;
    std::uint64_t steps;

    std::array<std::optional<ScaledEnclosure>, 126> values;
    P1Data(OwnedGeometry owner, QueryRequest original, LocalAmbientDescriptor layout,
           P1Accuracy request, P1Branch selected, P1EvaluationDiagnostic evaluation,
           P1WorkPlan work, std::uint64_t used) noexcept
        : geometry(std::move(owner)),
          query(original),
          ambient(layout),
          accuracy(request),
          branch(selected),
          diagnostic(evaluation),
          plan(work),
          steps(used) {}
};
struct P1Storage {
    std::optional<P1Data> data;
};
namespace p1 {
bool full_width_at_most(double radius_value, int output_scale, double tau,
                        const ExactInteger& longest, bool derivative, ExactProductScratch& product,
                        LocalRootScratch& scratch) {
    // Desired full width times ell or ell^2; square both sides at original SI exponents.
    auto& radius = scratch.integers[0];
    auto& tolerance = scratch.integers[1];
    auto& left = scratch.integers[2];
    auto& right = scratch.integers[3];
    if (radius_value == 0)
        return true;
    if (tau == 0)
        return false;
    local_decode(radius_value, radius);
    local_decode(tau, tolerance);
    // Remove exact q-grid powers before multiplying. Their thousands of zero
    // limbs are exponents, not arithmetic payload or a reason for heap scratch.
    const auto radius_power = boost::multiprecision::lsb(radius);
    const auto tau_power = boost::multiprecision::lsb(tolerance);
    const auto length_power = boost::multiprecision::lsb(longest);
    radius >>= radius_power;
    tolerance >>= tau_power;
    auto& length = scratch.integers[4];
    length = longest >> length_power;
    left = p1::exact_product(radius, radius, product);
    left = p1::exact_product(left, length, product);
    if (derivative)
        left = p1::exact_product(left, length, product);
    right = p1::exact_product(tolerance, tolerance, product);
    const int scale = derivative ? -8592 : -4296;
    const auto lp = static_cast<std::int64_t>(length_power);
    const int left_power =
        p1::exponent(2 * static_cast<std::int64_t>(output_scale) - 2146 + scale +
                     2 * static_cast<std::int64_t>(radius_power) + (derivative ? 2 : 1) * lp);
    const int right_power = p1::exponent(-2148 + 2 * static_cast<std::int64_t>(tau_power));
    return local_compare(left, left_power, right, right_power) <= 0;
}
double finite_upper(const Bound& bound) {
    if (bound.mantissa == 0)
        return 0;
    if (bound.mantissa < 0)
        throw ArithmeticFailure{core::FailureCode::PrecisionUnresolved};
    const int length = static_cast<int>(bits(bound.mantissa));
    const int shift = std::max({0, length - 53, -1074 - bound.power});
    if (shift >= length)
        return std::numeric_limits<double>::denorm_min();
    const Integer head = bound.mantissa >> static_cast<unsigned>(shift);
    double value = std::scalbn(head.convert_to<double>(), bound.power + shift);
    if ((head << static_cast<unsigned>(shift)) != bound.mantissa)
        value = std::nextafter(value, std::numeric_limits<double>::infinity());
    if (!std::isfinite(value))
        throw ArithmeticFailure{core::FailureCode::PrecisionUnresolved};
    return value;
}
P1QuadratureDifference qualify_far(std::span<const Interval, 126> low,
                                   std::span<const Interval, 126> high, const Interval& ell,
                                   int scale, Work& work) {
    P1QuadratureDifference differences{};
    const auto ell2 = multiply(ell, ell, work);
    for (std::size_t basis = 0; basis < 3; ++basis)
        for (std::size_t jet = 0; jet < 14; ++jet) {
            auto norm_difference = integer(0), norm_high = integer(0);
            for (std::size_t axis = 0; axis < 3; ++axis) {
                const auto index = 14 * (3 * basis + axis) + jet;
                const auto delta =
                    multiply(p1::scale(subtract(high[index], low[index], work), -scale),
                             jet == 0 ? ell : ell2, work);
                const auto c =
                    multiply(p1::scale(high[index], -scale), jet == 0 ? ell : ell2, work);
                norm_difference = add(norm_difference, multiply(delta, delta, work), work);
                norm_high = add(norm_high, multiply(c, c, work), work);
            }
            const auto delta = square_root(norm_difference, work);
            const auto c = square_root(norm_high, work);
            const auto absolute = word(1e-10), relative = multiply(word(1e-8), c, work);
            const auto threshold =
                compare(absolute.lower, relative.lower) > 0 ? absolute.lower : relative.lower;
            if (compare(delta.upper, threshold) > 0)
                throw ArithmeticFailure{core::FailureCode::PrecisionUnresolved};
            const double observed = finite_upper(delta.upper);
            if (jet == 0)
                differences.basis[basis] = observed;
            else
                differences.jacobian_columns[13 * basis + jet - 1] = observed;
        }
    return differences;
}
} // namespace p1
class P1Engine {
    static void point(const OwnedGeometry& geometry, PointId id, std::array<ExactInteger, 3>& out,
                      P1Scratch& scratch, p1::Work& work) {
        const auto coordinates = needed(geometry.point_words(id, Configuration::Current));
        for (std::size_t axis = 0; axis < 3; ++axis) {
            decode_words(coordinates[axis], out[axis], scratch.temporary[0], work);
            out[axis] = out[axis] < 0 ? -ExactInteger{(-out[axis]) << 1074}
                                      : ExactInteger{out[axis] << 1074};
        }
    }
    static void load(const OwnedGeometry& geometry, const QueryRequest& query,
                     const LocalAmbientDescriptor& ambient, P1Scratch& scratch, p1::Work& work) {
        for (std::size_t i = 0; i < 3; ++i)
            point(geometry, ambient.triangle_vertices[i], scratch.points[i], scratch, work);
        if (const auto* free = std::get_if<FreePoint>(&query.definition)) {
            for (std::size_t axis = 0; axis < 3; ++axis) {
                ScalarWords words = SingleWord{0};
                if (const auto* one = std::get_if<FreeSingle>(&free->words))
                    words = SingleWord{one->words[axis]};
                else {
                    const auto& pair = std::get<FreeDouble>(free->words);
                    words = DoubleWord{pair.high[axis], pair.low[axis]};
                }
                decode_words(words, scratch.points[3][axis], scratch.temporary[0], work);
                auto& coordinate = scratch.points[3][axis];
                coordinate = coordinate < 0 ? -ExactInteger{(-coordinate) << 1074}
                                            : ExactInteger{coordinate << 1074};
            }
        } else if (const auto* vertex = std::get_if<SourceVertex>(&query.definition))
            point(geometry, vertex->point, scratch.points[3], scratch, work);
        else {
            std::array<PointId, 3> vertices = ambient.triangle_vertices;
            std::size_t count = 0;
            scratch.weights.fill(0);
            scratch.weights[0] = ExactInteger{1} << 1074;
            if (const auto* edge = std::get_if<SourceEdgePoint>(&query.definition)) {
                const auto& data = *geometry.data_;
                const auto found = std::lower_bound(
                    data.edge_lookup.values().begin(), data.edge_lookup.values().end(),
                    edge->edge.value(),
                    [](const IdLookup& item, std::uint64_t id) { return item.id < id; });
                if (found == data.edge_lookup.values().end() || found->id != edge->edge.value())
                    throw p1::ArithmeticFailure{core::FailureCode::InvalidInput};
                const auto& source = data.edges.values()[static_cast<std::size_t>(found->ordinal)];
                for (std::size_t i = 0; i < 2; ++i)
                    vertices[i] = PointId{
                        data.points.values()[static_cast<std::size_t>(source.points[i])].id};
                decode_words(edge->t, scratch.weights[1], scratch.temporary[0], work);
                scratch.weights[0] -= scratch.weights[1];
                count = 2;
            } else if (const auto* face = std::get_if<SourceFacePoint>(&query.definition)) {
                vertices = needed(geometry.facet_points(face->facet));
                decode_words(face->u, scratch.weights[1], scratch.temporary[0], work);
                decode_words(face->v, scratch.weights[2], scratch.temporary[0], work);
                scratch.weights[0] -= scratch.weights[1] + scratch.weights[2];
                count = 3;
            }
            for (std::size_t axis = 0; axis < 3; ++axis) {
                scratch.points[3][axis] = 0;
                for (std::size_t i = 0; i < count; ++i) {
                    const auto coordinates =
                        needed(geometry.point_words(vertices[i], Configuration::Current));
                    decode_words(coordinates[axis], scratch.temporary[1], scratch.temporary[0],
                                 work);
                    work.step();
                    scratch.temporary[2] = p1::exact_product(scratch.temporary[1],
                                                             scratch.weights[i], scratch.product);
                    scratch.points[3][axis] += scratch.temporary[2];
                }
            }
        }
        decode_words(needed(geometry.core_radius()), scratch.temporary[3], scratch.temporary[0],
                     work);
    }
    static std::pair<p1::Options, int> prepare_input(P1Scratch& scratch, p1::Work& work) {
        p1::Options options;
        ExactInteger longest = 0;
        for (std::size_t edge = 0; edge < 3; ++edge) {
            ExactInteger length = 0, start = 0, end = 0;
            for (std::size_t axis = 0; axis < 3; ++axis) {
                work.step();
                const ExactInteger delta =
                    scratch.points[(edge + 1) % 3][axis] - scratch.points[edge][axis];
                const ExactInteger r = scratch.points[3][axis] - scratch.points[edge][axis];
                length += p1::exact_product(delta, delta, scratch.product);
                start += p1::exact_product(delta, r, scratch.product);
                end += p1::exact_product(delta, ExactInteger{r - delta}, scratch.product);
            }
            longest = std::max(longest, length);
            options.start_sign[edge] = sign(start);
            options.end_sign[edge] = sign(end);
        }
        const int k = static_cast<int>(boost::multiprecision::msb(longest)) - 4296;
        const int scale = k >= 0 ? k / 2 : -((-k + 1) / 2);
        std::array<ExactInteger, 3> e, f, r;
        for (std::size_t axis = 0; axis < 3; ++axis) {
            e[axis] = scratch.points[1][axis] - scratch.points[0][axis];
            f[axis] = scratch.points[2][axis] - scratch.points[0][axis];
            r[axis] = scratch.points[3][axis] - scratch.points[0][axis];
        }
        ExactInteger plane = 0;
        for (std::size_t axis = 0; axis < 3; ++axis) {
            const auto next = (axis + 1) % 3, last = (axis + 2) % 3;
            const ExactInteger cross = p1::exact_product(e[next], f[last], scratch.product) -
                                       p1::exact_product(e[last], f[next], scratch.product);
            plane += p1::exact_product(cross, r[axis], scratch.product);
        }
        options.exterior = plane == 0 && scratch.temporary[3] == 0;
        for (std::size_t axis = 0; axis < 3; ++axis) {
            scratch.input[axis] = p1::exact(r[axis], -2148 - scale);
            for (std::size_t vertex = 0; vertex < 3; ++vertex)
                scratch.input[3 + 3 * vertex + axis] = p1::exact(
                    scratch.points[vertex][axis] - scratch.points[0][axis], -2148 - scale);
        }
        scratch.input[12] = p1::exact(scratch.temporary[3], -1074 - scale);
        // Reuse exact longest for final full-width tests after all query scratch is done.
        scratch.temporary[4] = longest;
        return {options, scale};
    }
    static double upward(const ExactInteger& value, int power) {
        if (value == 0)
            return 0;
        const int length = static_cast<int>(boost::multiprecision::msb(value)) + 1;
        const int shift = std::max({0, length - 53, -1074 - power});
        if (shift >= length)
            return std::numeric_limits<double>::denorm_min();
        const ExactInteger head = value >> static_cast<unsigned>(shift);
        double result = std::scalbn(head.convert_to<double>(), power + shift);
        if ((head << static_cast<unsigned>(shift)) != value)
            result = std::nextafter(result, std::numeric_limits<double>::infinity());
        if (!std::isfinite(result))
            throw p1::ArithmeticFailure{core::FailureCode::PrecisionUnresolved};
        return result;
    }

public:
    static ScaledEnclosure enclosure(const p1::Interval& interval, core::Precision precision,
                                     LocalRootScratch& scratch) {
        const auto top = std::max(
            static_cast<std::int64_t>(p1::bits(interval.lower.mantissa)) + interval.lower.power,
            static_cast<std::int64_t>(p1::bits(interval.upper.mantissa)) + interval.upper.power);
        const int common = p1::exponent(top - 258);
        const auto low = p1::at_power(interval.lower, common, false),
                   high = p1::at_power(interval.upper, common, true);
        const p1::Integer middle = (low + high) / 2;
        const auto center =
            LocalBounds::dyadic(middle.convert_to<ExactInteger>(), common, precision, scratch);
        if (!center)
            throw p1::ArithmeticFailure{center.error().code};
        const p1::Integer width = std::max(middle - low, high - middle);
        const double additional =
            upward(width.convert_to<ExactInteger>(), common - center->exponent());
        const double radius =
            std::nextafter(center->radius() + additional, std::numeric_limits<double>::infinity());
        const double exact_radius = (center->radius() == 0 && additional == 0) ? 0 : radius;
        if (!std::isfinite(exact_radius))
            throw p1::ArithmeticFailure{core::FailureCode::PrecisionUnresolved};
        return ScaledEnclosure{center->center(), exact_radius, center->exponent()};
    }

private:
    static bool width(const ScaledEnclosure& value, double tau, const ExactInteger& longest,
                      bool derivative, P1Scratch& scratch) {
        return p1::full_width_at_most(value.radius(), value.exponent(), tau, longest, derivative,
                                      scratch.product, scratch.output);
    }

public:
    static AeroResult<P1Basis> evaluate(const LocalGeometryPacket&, P1Accuracy, P1Limits,
                                        P1Workspace&, core::ResourcePlan&, const std::stop_token&,
                                        const std::array<bool, 3>&) noexcept;
};
std::optional<core::FailureCode> P1Accuracy::validation_error() const noexcept {
    if (!std::isfinite(tau_basis) || !std::isfinite(tau_jacobian))
        return core::FailureCode::Nonfinite;
    if (tau_basis < 0 || tau_jacobian < 0)
        return core::FailureCode::InvalidInput;
    if (representation != core::Precision::Binary64 && representation != core::Precision::TwoWord)
        return core::FailureCode::PrecisionUnsupported;
    return std::nullopt;
}

AeroResult<P1Basis> P1Engine::evaluate(const LocalGeometryPacket& packet, P1Accuracy accuracy,
                                       P1Limits limits, P1Workspace& workspace,
                                       core::ResourcePlan& plan, const std::stop_token& stop,
                                       const std::array<bool, 3>& boundary) noexcept {
    const auto context = packet.context();
    if (workspace.entered_.test_and_set())
        return std::unexpected(failed(core::FailureCode::InvalidInput, context));
    const struct Release {
        std::atomic_flag& flag;
        ~Release() {
            flag.clear();
        }
    } release{workspace.entered_};
    try {
        if (!supported_environment())
            return std::unexpected(failed(core::FailureCode::PrecisionUnsupported, context));
        if (!packet.active())
            return std::unexpected(failed(core::FailureCode::InvalidInput, context));
        if (const auto error = accuracy.validation_error())
            return std::unexpected(failed(*error, context));
        if (stop.stop_requested())
            return std::unexpected(failed(core::FailureCode::Cancelled, context));
        const auto mode = limits.external_memory.selection();
        if (std::holds_alternative<numerics::linear::RequireProvenBound>(mode))
            return std::unexpected(failed(core::FailureCode::Unsupported, context));
        if (std::get<numerics::linear::EstimateAdmission>(mode).bytes <
            predicate_memory_plan().peak_payload_estimate_bytes())
            return std::unexpected(failed(core::FailureCode::ResourceExceeded, context));
        // The immutable expression plan belongs to this prepared workspace.
        const Shape shape = workspace.scratch_.values()[0].shape;
        if (limits.max_steps < shape.work || plan.remaining_bytes() < sizeof(P1Storage))
            return std::unexpected(failed(core::FailureCode::ResourceExceeded, context));
        const auto& geometry = needed(packet.geometry()).get();
        if (geometry.context() != context ||
            needed(geometry.identity()) != needed(packet.snapshot_identity()))
            return std::unexpected(failed(core::FailureCode::StaleGeneration, context));
        const auto query = needed(packet.query());
        const auto ambient = needed(packet.ambient());
        const auto far = needed(packet.far());
        auto reserved = plan.allocate<P1Storage>(1, operation);
        if (!reserved)
            return std::unexpected(OwnedGeometryFailure{reserved.error(), {}, context});
        // User allocation hooks may change the thread environment during reservation.
        if (!supported_environment())
            return std::unexpected(failed(core::FailureCode::PrecisionUnsupported, context));
        auto& scratch = workspace.scratch_.values()[0];
        p1::Work work{0, limits.max_steps, stop};
        load(geometry, query, ambient, scratch, work);
        auto [options, scale] = prepare_input(scratch, work);
        options.boundary = boundary;
        scratch.pi = p1::pi_constant(work);
        scratch.ln2 = p1::logarithm_series(p1::divide(p1::integer(1), p1::integer(3), work), work);
        p1::Program program{scratch.program.values()};
        const auto output =
            far == FarRelation::Far ? p1::far_program(program) : p1::near_program(program, options);
        if (program.count > scratch.jets.size())
            return std::unexpected(failed(core::FailureCode::ResourceExceeded, context));
        const auto evaluate = [&]() {
            for (std::size_t i = 0; i < program.count; ++i)
                p1::evaluate_node(scratch.program.values()[i], scratch.jets.values(), i,
                                  scratch.input, -scale, scratch.pi, scratch.ln2, work);
        };
        P1EvaluationDiagnostic diagnostic = P1AnalyticEvaluation{};
        const std::size_t base = program.count;
        if (far == FarRelation::Far) {
            p1::build_rule(p1::gauss_8, scratch.low, work);
            p1::build_rule(p1::gauss_12, scratch.high, work);
            for (std::size_t i = 0; i < 18; ++i)
                scratch.jets.values()[base + i].components.fill(p1::integer(0));
            for (std::size_t order = 0; order < 2; ++order) {
                const auto& quadrature = order == 0 ? scratch.low : scratch.high;
                for (std::size_t i = 0; i < quadrature.count; ++i)
                    for (std::size_t j = 0; j < quadrature.count; ++j) {
                        scratch.input[13] = quadrature.nodes[i];
                        scratch.input[14] = quadrature.nodes[j];
                        scratch.input[15] =
                            p1::multiply(quadrature.weights[i], quadrature.weights[j], work);
                        evaluate();
                        for (std::size_t component = 0; component < 9; ++component)
                            for (std::size_t jet = 0; jet < 14; ++jet) {
                                auto& accumulator =
                                    scratch.jets.values()[base + 9 * order + component]
                                        .components[jet];
                                accumulator = p1::add(
                                    accumulator,
                                    scratch.jets.values()[output.nodes[component]].components[jet],
                                    work);
                            }
                    }
            }
            for (std::size_t component = 0; component < 9; ++component)
                for (std::size_t jet = 0; jet < 14; ++jet) {
                    scratch.low_values[14 * component + jet] =
                        scratch.jets.values()[base + component].components[jet];
                    scratch.high_values[14 * component + jet] =
                        scratch.jets.values()[base + 9 + component].components[jet];
                }
            const auto ell = p1::square_root(p1::exact(scratch.temporary[4], -4296), work);
            diagnostic = p1::qualify_far(scratch.low_values, scratch.high_values, ell, scale, work);
        } else
            evaluate();
        const auto selected =
            far == FarRelation::Far
                ? P1Branch::FarGauss12
                : (options.exterior ? P1Branch::IdealCoplanarExterior : P1Branch::Analytic);
        const P1WorkPlan admitted{shape.nodes, shape.work,
                                  p1::checked_sum(shape_payload(shape), sizeof(P1Storage))};
        auto& result = reserved->values()[0].data.emplace(OwnedGeometry{geometry.data_, context},
                                                          query, ambient, accuracy, selected,
                                                          diagnostic, admitted, 0);
        for (std::size_t component = 0; component < 9; ++component) {
            const auto& evaluated =
                scratch.jets.values()[far == FarRelation::Far ? base + 9 + component
                                                              : output.nodes[component]];
            for (std::size_t jet = 0; jet < 14; ++jet) {
                work.step();
                auto value = enclosure(p1::scale(evaluated.components[jet], -scale),
                                       accuracy.representation, scratch.output);
                if (!width(value, jet == 0 ? accuracy.tau_basis : accuracy.tau_jacobian,
                           scratch.temporary[4], jet != 0, scratch))
                    throw p1::ArithmeticFailure{core::FailureCode::PrecisionUnresolved};
                const std::size_t index = jet == 0 ? component : 9 + 13 * component + jet - 1;
                result.values[index] = value;
            }
        }
        result.steps = work.used;
        return P1Basis{std::move(*reserved), context};
    } catch (const OwnedGeometryFailure& failure) {
        return std::unexpected(failure);
    } catch (const p1::ArithmeticFailure& failure) {
        return std::unexpected(failed(failure.code, context));
    } catch (const std::bad_alloc&) {
        return std::unexpected(failed(core::FailureCode::ResourceExceeded, context));
    } catch (...) {
        return std::unexpected(failed(core::FailureCode::PrecisionUnresolved, context));
    }
}
P1Workspace::P1Workspace(Key, core::PlannedBuffer<P1Scratch> scratch) noexcept
    : scratch_(std::move(scratch)) {}
P1Workspace::~P1Workspace() = default;
core::Result<std::unique_ptr<P1Workspace>>
P1Workspace::create(core::ResourcePlan& plan, core::FailureDetails details) noexcept {
    try {
        const auto shape = maximum_shape();
        const auto bytes = shape_payload(shape);
        if (plan.remaining_bytes() < bytes)
            return std::unexpected(
                core::Failure{core::FailureCode::ResourceExceeded, operation, details});
        auto scratch = plan.allocate<P1Scratch>(1, operation);
        if (!scratch)
            return std::unexpected(with_details(scratch.error(), details));
        auto nodes = plan.allocate<p1::Node>(static_cast<core::Index>(shape.nodes), operation);
        if (!nodes)
            return std::unexpected(with_details(nodes.error(), details));
        auto jets = plan.allocate<p1::Jet>(static_cast<core::Index>(shape.nodes), operation);
        if (!jets)
            return std::unexpected(with_details(jets.error(), details));
        scratch->values()[0].shape = shape;
        scratch->values()[0].program = std::move(*nodes);
        scratch->values()[0].jets = std::move(*jets);
        return std::make_unique<P1Workspace>(Key{}, std::move(*scratch));
    } catch (const p1::ArithmeticFailure& failure) {
        return std::unexpected(core::Failure{failure.code, operation, details});
    } catch (const std::bad_alloc&) {
        return std::unexpected(
            core::Failure{core::FailureCode::ResourceExceeded, operation, details});
    } catch (...) {
        return std::unexpected(
            core::Failure{core::FailureCode::PrecisionUnresolved, operation, details});
    }
}
std::uint64_t P1Workspace::payload_bytes() noexcept {
    try {
        return shape_payload(maximum_shape());
    } catch (...) {
        return std::numeric_limits<std::uint64_t>::max();
    }
}
P1Basis::P1Basis(core::PlannedBuffer<P1Storage> storage, GeometryContext context) noexcept
    : storage_(std::move(storage)), context_(context) {}
P1Basis::P1Basis(P1Basis&&) noexcept = default;
P1Basis& P1Basis::operator=(P1Basis&&) noexcept = default;
P1Basis::~P1Basis() = default;
const P1Data* P1Basis::data() const noexcept {
    if (storage_.size() != 1)
        return nullptr;
    const auto& owned = storage_.values()[0].data;
    if (!owned)
        return nullptr;
    return std::addressof(*owned);
}
bool P1Basis::active() const noexcept {
    return data() != nullptr;
}
AeroResult<std::reference_wrapper<const ScaledEnclosure>>
P1Basis::basis(std::size_t k, std::size_t xyz) const& noexcept {
    const auto* owned = data();
    if (owned == nullptr || k >= 3 || xyz >= 3)
        return std::unexpected(failed(core::FailureCode::InvalidInput, context_));
    const auto& value = owned->values[3 * k + xyz];
    if (!value)
        return std::unexpected(failed(core::FailureCode::InvalidInput, context_));
    return std::cref(*value);
}
AeroResult<std::reference_wrapper<const ScaledEnclosure>>
P1Basis::derivative(std::size_t k, std::size_t xyz, std::size_t variable) const& noexcept {
    const auto* owned = data();
    if (owned == nullptr || k >= 3 || xyz >= 3 || variable >= 13)
        return std::unexpected(failed(core::FailureCode::InvalidInput, context_));
    const auto& value = owned->values[9 + 13 * (3 * k + xyz) + variable];
    if (!value)
        return std::unexpected(failed(core::FailureCode::InvalidInput, context_));
    return std::cref(*value);
}
AeroResult<std::reference_wrapper<const OwnedGeometry>> P1Basis::geometry() const& noexcept {
    const auto* owned = data();
    if (owned == nullptr)
        return std::unexpected(failed(core::FailureCode::InvalidInput, context_));
    return std::cref(owned->geometry);
}
AeroResult<QueryRequest> P1Basis::query() const noexcept {
    const auto* owned = data();
    if (owned == nullptr)
        return std::unexpected(failed(core::FailureCode::InvalidInput, context_));
    return owned->query;
}
AeroResult<LocalAmbientDescriptor> P1Basis::ambient() const noexcept {
    const auto* owned = data();
    if (owned == nullptr)
        return std::unexpected(failed(core::FailureCode::InvalidInput, context_));
    return owned->ambient;
}
AeroResult<P1Branch> P1Basis::branch() const noexcept {
    const auto* owned = data();
    if (owned == nullptr)
        return std::unexpected(failed(core::FailureCode::InvalidInput, context_));
    return owned->branch;
}
AeroResult<P1EvaluationDiagnostic> P1Basis::diagnostic() const noexcept {
    const auto* owned = data();
    if (owned == nullptr)
        return std::unexpected(failed(core::FailureCode::InvalidInput, context_));
    return owned->diagnostic;
}
AeroResult<P1Accuracy> P1Basis::accuracy() const noexcept {
    const auto* owned = data();
    if (owned == nullptr)
        return std::unexpected(failed(core::FailureCode::InvalidInput, context_));
    return owned->accuracy;
}
AeroResult<std::uint64_t> P1Basis::steps() const noexcept {
    const auto* owned = data();
    if (owned == nullptr)
        return std::unexpected(failed(core::FailureCode::InvalidInput, context_));
    return owned->steps;
}
AeroResult<P1WorkPlan> P1Basis::work_plan() const noexcept {
    const auto* owned = data();
    if (owned == nullptr)
        return std::unexpected(failed(core::FailureCode::InvalidInput, context_));
    return owned->plan;
}
AeroResult<std::array<std::byte, 32>> P1Basis::snapshot_identity() const noexcept {
    const auto* owned = data();
    if (owned == nullptr)
        return std::unexpected(failed(core::FailureCode::InvalidInput, context_));
    return owned->geometry.identity();
}
std::uint64_t P1Basis::payload_bytes() noexcept {
    return sizeof(P1Storage);
}
AeroResult<P1Basis> evaluate_p1_basis(const LocalGeometryPacket& packet, P1Accuracy accuracy,
                                      P1Limits limits, P1Workspace& workspace,
                                      core::ResourcePlan& plan,
                                      const std::stop_token& stop) noexcept {
    return P1Engine::evaluate(packet, accuracy, limits, workspace, plan, stop, {true, true, true});
}

bool p1_environment_supported() noexcept {
    return supported_environment();
}

p1::Interval decode_p1_word(double value) {
    return word(value);
}

ScaledEnclosure enclose_p1_value(const p1::Interval& value, core::Precision precision,
                                 LocalRootScratch& scratch) {
    return P1Engine::enclosure(value, precision, scratch);
}

AeroResult<P1Basis> evaluate_p1_partition(const LocalGeometryPacket& packet,
                                          const std::array<bool, 3>& boundary, P1Accuracy accuracy,
                                          P1Limits limits, P1Workspace& workspace,
                                          core::ResourcePlan& plan,
                                          const std::stop_token& stop) noexcept {
    return P1Engine::evaluate(packet, accuracy, limits, workspace, plan, stop, boundary);
}
} // namespace uvlas::aero::detail
