#include "pg_graph_detail.hpp"

#include <boost/hash2/sha2.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <charconv>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <string_view>
#include <system_error>
#include <type_traits>
#include <utility>

namespace uvlas::structure {
namespace {
constexpr auto operation = core::Operation::PrepareModel;
using detail::IntegerBackend;
using detail::PgCellPointData;
using detail::PgCellPolynomialsData;
using detail::PgGraphData;
using detail::PgGraphFactory;
using detail::PgGraphValueData;
using detail::PgRectangleActivityData;
struct Stop {
    core::Failure failure;
};
[[noreturn]] void fail(core::FailureCode code) {
    throw Stop{core::Failure{code, operation}};
}
std::uint64_t add(std::uint64_t a, std::uint64_t b) {
    if (b > std::numeric_limits<std::uint64_t>::max() - a)
        fail(core::FailureCode::ResourceExceeded);
    return a + b;
}
std::uint64_t mul(std::uint64_t a, std::uint64_t b) {
    if (a != 0 && b > std::numeric_limits<std::uint64_t>::max() / a)
        fail(core::FailureCode::ResourceExceeded);
    return a * b;
}
core::Index index(std::uint64_t v) {
    if (v > static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max()))
        fail(core::FailureCode::ResourceExceeded);
    return static_cast<core::Index>(v);
}
template <class T>
T take(core::Result<T> result) {
    if (!result)
        throw Stop{result.error()};
    return std::move(*result);
}
void cancelled(const PgGraphLimits& limits) {
    if (limits.cancellation.stop_requested())
        fail(core::FailureCode::Cancelled);
}
void validate(const PgGraphLimits& limits) {
    const auto& a = limits.algebra;
    if (a.external_memory.valueless_by_exception() ||
        (a.timing != ExactTimingMode::Disabled && a.timing != ExactTimingMode::Enabled))
        fail(core::FailureCode::InvalidInput);
    for (const auto value :
         {a.maximum_input_characters, a.maximum_terms, a.maximum_raw_terms, a.maximum_axis_degree,
          a.maximum_total_degree, a.maximum_value_bits, a.maximum_work_bits, limits.maximum_cells})
        (void)index(value);
}
ExactLimits inner(const PgGraphLimits& limits) {
    auto result = limits.algebra;
    result.timing = ExactTimingMode::Disabled;
    return result;
}
PgGraphLimits nested(const PgGraphLimits& limits) {
    return {inner(limits), limits.maximum_cells, limits.cancellation};
}
std::uint64_t value_bits(const ExactRational& value) noexcept {
    return std::max(value.numerator_bits(), value.denominator_bits());
}
void check_value(const ExactRational& value, const PgGraphLimits& limits) {
    if (value_bits(value) > limits.algebra.maximum_value_bits)
        fail(core::FailureCode::ResourceExceeded);
}
std::uint64_t retained(const ExactRational& value) {
    const auto& b = detail::ExactFactory::value(value).backend();
    const auto capacity = [](const IntegerBackend& integer) {
        return integer.capacity() > IntegerBackend::internal_limb_count ? integer.capacity() : 0;
    };
    return mul(add(capacity(b.num()), capacity(b.denom())),
               sizeof(boost::multiprecision::limb_type));
}
std::uint64_t retained(const PgGraphData& data) {
    auto total = add(add(retained(data.x0), retained(data.width)), retained(data.skew));
    for (const auto& term : data.height.terms())
        total = add(total, retained(term.coefficient));
    return total;
}
std::uint64_t maximum_bits(const ExactChartPolynomial& f, const PgGraphLimits& limits) {
    std::uint64_t result = 1;
    if (f.terms().size() > limits.algebra.maximum_terms)
        fail(core::FailureCode::ResourceExceeded);
    for (const auto& term : f.terms()) {
        if (term.powers.first < 0 || term.powers.second < 0)
            fail(core::FailureCode::InvalidInput);
        const auto a = static_cast<std::uint64_t>(term.powers.first),
                   b = static_cast<std::uint64_t>(term.powers.second);
        if (a > limits.algebra.maximum_axis_degree || b > limits.algebra.maximum_axis_degree ||
            add(a, b) > limits.algebra.maximum_total_degree)
            fail(core::FailureCode::ResourceExceeded);
        check_value(term.coefficient, limits);
        result = std::max(result, value_bits(term.coefficient));
    }
    return result;
}
void definition(const PgGraphRequest& request, const PgGraphLimits& limits) {
    if (request.n != 1 && request.n != 2 && request.n != 4)
        fail(core::FailureCode::InvalidInput);
    if (request.convention != PgGraphConvention::Version1 || request.width.sign() <= 0)
        fail(core::FailureCode::InvalidInput);
    const auto degree = request.height.degree();
    if (degree && *degree > 4)
        fail(core::FailureCode::InvalidInput);
    check_value(request.x0, limits);
    check_value(request.width, limits);
    check_value(request.skew, limits);
    (void)maximum_bits(request.height, limits);
    if (mul(static_cast<std::uint64_t>(2 * request.n), static_cast<std::uint64_t>(2 * request.n)) >
        limits.maximum_cells)
        fail(core::FailureCode::ResourceExceeded);
}
const PgGraphData& graph_data(const PgGraph& graph, const PgGraphLimits& limits) {
    const auto* data = PgGraphFactory::data(graph);
    if (data == nullptr || !data->identity)
        fail(core::FailureCode::InvalidInput);
    definition({data->x0, data->width, data->skew, data->height, data->n, data->convention},
               limits);
    return *data;
}
const PgGraphContentIdentity& identity(const PgGraphData& data) {
    const auto& value = data.identity;
    if (!value)
        fail(core::FailureCode::InvalidInput);
    return *value;
}
std::uint64_t graph_bits(const PgGraphData& d, const PgGraphLimits& l) {
    return std::max(
        {value_bits(d.x0), value_bits(d.width), value_bits(d.skew), maximum_bits(d.height, l)});
}
bool contents_equal(const PgGraphData& a, const PgGraphData& b) noexcept {
    if (a.n != b.n || a.convention != b.convention || !(a.x0 == b.x0) || !(a.width == b.width) ||
        !(a.skew == b.skew))
        return false;
    const auto x = a.height.terms(), y = b.height.terms();
    if (x.size() != y.size())
        return false;
    for (std::size_t i = 0; i < x.size(); ++i)
        if (x[i].powers.first != y[i].powers.first || x[i].powers.second != y[i].powers.second ||
            !(x[i].coefficient == y[i].coefficient))
            return false;
    const auto ca = a.cells.values(), cb = b.cells.values();
    return ca.size() == cb.size() && std::equal(ca.begin(), ca.end(), cb.begin());
}
void check_cell(const PgGraphData& d, PgCellIndex cell) {
    if (cell.first < 0 || cell.second < 0 || cell.first >= 2 * d.n || cell.second >= 2 * d.n)
        fail(core::FailureCode::InvalidInput);
}
void check_side(PgSide side) {
    if (side != PgSide::Interior && side != PgSide::Lower && side != PgSide::Upper)
        fail(core::FailureCode::InvalidInput);
}
void check_orders(PgDerivativeOrders orders, PgDerivativeFrame frame) {
    if (orders.first < 0 || orders.second < 0 ||
        (frame != PgDerivativeFrame::ChartST && frame != PgDerivativeFrame::CellUV))
        fail(core::FailureCode::InvalidInput);
}
bool structural_zero(const PgGraphData& d, PgDerivativeOrders orders) noexcept {
    // Test before adding orders or forming a scale: even maximal Index orders are valid.
    if (orders.first > 4 || orders.second > 4 || orders.first + orders.second > 4)
        return true;
    const auto total = orders.first + orders.second;
    if (total <= 1)
        return false; // x_s=w>0 and y_t=1.
    if (d.height.terms().empty())
        return true;
    return std::none_of(
        d.height.terms().begin(), d.height.terms().end(), [&](const ExactTerm& term) {
            return term.powers.first >= orders.first && term.powers.second >= orders.second;
        });
}
std::uint64_t poly_payload(std::uint64_t count) {
    const auto item = add(add(detail::exact_polynomial_raw_term_bytes(), sizeof(ExactTerm)),
                          sizeof(detail::RationalStorage));
    return add(mul(count, item), sizeof(detail::PolynomialStorage));
}
std::uint64_t prepare_payload(const PgGraphData& d) {
    return add(add(sizeof(PgGraphData), mul(d.cells.size(), sizeof(PgCellIndex))),
               add(mul(3, sizeof(detail::RationalStorage)), poly_payload(d.height.terms().size())));
}
struct Bounds {
    Bounds(std::uint64_t work_bound, std::uint64_t raw_bound, std::uint64_t held_bound,
           std::uint64_t payload_bound, std::uint64_t input_payload, bool copying = false) noexcept
        : work(work_bound),
          raw(raw_bound),
          held(held_bound),
          payload(payload_bound),
          borrowed(input_payload),
          copy_only(copying) {}
    std::uint64_t work, raw, held, payload, borrowed;
    bool copy_only;
};
PgGraphCallInfo admit(const Bounds& bounds, const PgGraphLimits& limits, core::ResourcePlan& plan) {
    if (bounds.raw > limits.algebra.maximum_raw_terms)
        fail(core::FailureCode::ResourceExceeded);
    if (std::holds_alternative<ExactRequireProvenBound>(limits.algebra.external_memory))
        fail(core::FailureCode::Unsupported);
    constexpr auto width = IntegerBackend::limb_bits;
    const auto limbs =
        bounds.work == 0 ? 0 : add(add(bounds.work / width, bounds.work % width != 0 ? 1 : 0), 2);
    const auto guarded = mul(limbs, width);
    if (guarded > limits.algebra.maximum_work_bits)
        fail(core::FailureCode::ResourceExceeded);
    const auto capacity = limbs <= IntegerBackend::internal_limb_count ? 0 : mul(4, limbs);
    const auto resize = capacity == 0 ? 0 : limbs;
    const auto gcd = limbs == 0 ? 0 : add(mul(6, add(limbs, 1)), mul(5, capacity));
    const auto divide = mul(2, capacity), karatsuba = limbs >= 60 ? mul(5, limbs) : 0;
    const auto slots = bounds.work == 0   ? 0
                       : bounds.copy_only ? 2
                                          : add(mul(2, bounds.held), add(16, mul(6, bounds.raw)));
    const auto external = add(
        bounds.borrowed, mul(add(add(mul(slots, capacity),
                                     (bounds.copy_only ? 0 : std::max({gcd, divide, karatsuba}))),
                                 resize),
                             sizeof(boost::multiprecision::limb_type)));
    const auto* permission = std::get_if<ExactEstimateAdmission>(&limits.algebra.external_memory);
    if (permission == nullptr)
        fail(core::FailureCode::InvalidInput);
    if (external > permission->bytes || bounds.payload > plan.remaining_bytes())
        fail(core::FailureCode::ResourceExceeded);
    return {{external, bounds.work, guarded, bounds.raw, std::nullopt}, bounds.payload, 0};
}
Bounds prepare_bounds(const PgGraphRequest& r, const PgGraphLimits& l) {
    const auto m = r.height.terms().size();
    const auto f = maximum_bits(r.height, l);
    auto borrowed = add(add(retained(r.x0), retained(r.width)), retained(r.skew));
    for (const auto& term : r.height.terms())
        borrowed = add(borrowed, retained(term.coefficient));
    const auto count =
        mul(static_cast<std::uint64_t>(2 * r.n), static_cast<std::uint64_t>(2 * r.n));
    return {
        std::max({value_bits(r.x0), value_bits(r.width), value_bits(r.skew), mul(m, add(f, 1))}), m,
        18,
        add(add(sizeof(PgGraphData), mul(count, sizeof(PgCellIndex))),
            add(mul(3, sizeof(detail::RationalStorage)), poly_payload(m))),
        borrowed};
}
void native_chars(std::uint64_t value, const PgGraphLimits& limits) {
    std::uint64_t count = 1;
    while (value >= 10) {
        value /= 10;
        ++count;
    }
    if (count > limits.algebra.maximum_input_characters)
        fail(core::FailureCode::ResourceExceeded);
}
Bounds point_bounds(const PgGraphData& d, const ExactRational& a, const ExactRational& b,
                    const PgGraphLimits& l) {
    native_chars(static_cast<std::uint64_t>(4 * d.n), l);
    const auto base = prepare_bounds({d.x0, d.width, d.skew, d.height, d.n, d.convention}, l);
    return {std::max({base.work, add(std::max(value_bits(a), value_bits(b)), 7), std::uint64_t{9}}),
            base.raw, 27,
            add(add(prepare_payload(d), sizeof(PgCellPointData)),
                mul(9, sizeof(detail::RationalStorage))),
            add(retained(d), add(retained(a), retained(b)))};
}
std::uint64_t substitution_raw(const ExactChartPolynomial& f) {
    std::uint64_t count = 0;
    for (const auto& term : f.terms())
        count = add(count, std::uint64_t{1}
                               << static_cast<unsigned>(term.powers.first + term.powers.second));
    return count;
}
Bounds cell_bounds(const PgGraphData& d, const PgGraphLimits& l) {
    if (l.algebra.maximum_axis_degree < 1 || l.algebra.maximum_total_degree < 1 ||
        l.algebra.maximum_terms < 2)
        fail(core::FailureCode::ResourceExceeded);
    native_chars(static_cast<std::uint64_t>(4 * d.n), l);
    const auto f = maximum_bits(d.height, l), v = graph_bits(d, l), u = substitution_raw(d.height);
    const auto sub =
        add(mul(3 * u - d.height.terms().size(), detail::exact_polynomial_raw_term_bytes()),
            add(mul(15, add(sizeof(ExactTerm), sizeof(detail::RationalStorage))),
                sizeof(detail::PolynomialStorage)));
    return {std::max({mul(u, add(f, 21)), mul(3, add(mul(3, v), 12)), std::uint64_t{12}}),
            std::max(u, std::uint64_t{3}), 35,
            add(add(sizeof(PgCellPolynomialsData), mul(11, sizeof(detail::RationalStorage))),
                add(add(mul(3, poly_payload(2)), poly_payload(3)), sub)),
            retained(d)};
}
Bounds evaluation_bounds(const PgGraphData& d, const PgCellPointData& p, PgDerivativeOrders orders,
                         PgDerivativeFrame frame, const PgGraphLimits& l) {
    const auto& source = p.source;
    if (!source)
        fail(core::FailureCode::InvalidInput);
    const auto& snapshot = graph_data(*source, l);
    auto borrowed = add(retained(d), retained(snapshot));
    borrowed =
        add(borrowed, add(add(retained(p.s), retained(p.t)), add(retained(p.u), retained(p.v))));
    if (structural_zero(d, orders))
        return {0, 0, 0, sizeof(PgGraphValueData), borrowed};
    const auto total = orders.first + orders.second;
    if (frame == PgDerivativeFrame::CellUV && total != 0) {
        std::uint64_t denominator = 1;
        for (core::Index k = 0; k < total; ++k)
            denominator = mul(denominator, static_cast<std::uint64_t>(4 * d.n));
        native_chars(denominator, l);
    } else if (orders.first == 0 && orders.second == 1)
        native_chars(1, l);
    const auto f = maximum_bits(d.height, l), v = graph_bits(d, l), m = d.height.terms().size(),
               point = std::max(value_bits(p.s), value_bits(p.t));
    const auto work =
        std::max({mul(m, add(f, 13)), add(mul(m, add(add(f, 12), add(mul(4, point), 1))), 17),
                  add(add(mul(3, v), mul(2, point)), 20), std::uint64_t{21}});
    return {work, m, 42,
            add(add(sizeof(PgGraphValueData), mul(8, sizeof(detail::RationalStorage))),
                mul(2, poly_payload(m))),
            borrowed};
}
PgGraphCallInfo completed(PgGraphCallInfo info, std::uint64_t before,
                          const core::ResourcePlan& plan) {
    info.allocated_owned_payload_bytes = plan.allocated_bytes() - before;
    // A layout/count regression must never be hidden behind a successful forecast.
    if (info.allocated_owned_payload_bytes > info.planned_owned_payload_upper_bytes)
        fail(core::FailureCode::ResourceExceeded);
    return info;
}
template <class F>
auto entry(const PgGraphLimits& limits, F&& callable) noexcept -> std::invoke_result_t<F> {
    using Result = std::invoke_result_t<F>;
    core::StageRecorder recorder{limits.algebra.timing == ExactTimingMode::Enabled, operation};
    auto result =
        core::guard_entry(operation, core::FailureCode::InvalidInput, {}, [&]() -> Result {
            const core::StageScope stage{recorder, core::MeasurementStage::Preparation};
            try {
                validate(limits);
                return callable();
            } catch (const Stop& stopped) {
                return std::unexpected(stopped.failure);
            }
        });
    const auto timing = recorder.snapshot();
    if (!timing)
        return std::unexpected(timing.error());
    if (result && recorder.enabled()) {
        const auto sample =
            (*timing)[static_cast<std::size_t>(core::MeasurementStage::Preparation)];
        if constexpr (std::is_same_v<typename Result::value_type, ExactRational>) {
            detail::ExactFactory::timing(*result, sample);
        } else {
            auto info = result->preparation();
            info.arithmetic.timing = sample;
            PgGraphFactory::info(*result, info);
        }
    }
    return result;
}
// Minimal big-endian magnitude streaming; no text or dynamic canonical byte buffer.
void append_u64(boost::hash2::sha2_256& hash, std::uint64_t value) noexcept {
    std::array<unsigned char, 8> bytes{};
    for (std::size_t i = 0; i < bytes.size(); ++i) {
        bytes[bytes.size() - 1 - i] = static_cast<unsigned char>(value & 0xffU);
        value >>= 8U;
    }
    hash.update(bytes.data(), bytes.size());
}
void append_text(boost::hash2::sha2_256& hash, std::string_view text) noexcept {
    append_u64(hash, text.size());
    if (!text.empty())
        hash.update(text.data(), text.size());
}
std::uint64_t backend_bits(const IntegerBackend& value) noexcept {
    if (value.size() == 1 && value.limbs()[0] == 0)
        return 0;
    return static_cast<std::uint64_t>(
        (value.size() - 1) * IntegerBackend::limb_bits + IntegerBackend::limb_bits -
        static_cast<std::size_t>(std::countl_zero(value.limbs()[value.size() - 1])));
}
void append_integer(boost::hash2::sha2_256& hash, const IntegerBackend& value) {
    const auto bytes = add(backend_bits(value), 7) / 8;
    append_u64(hash, bytes);
    for (auto i = bytes; i != 0; --i) {
        const auto bit = (i - 1) * 8;
        const auto byte = static_cast<unsigned char>(
            (value.limbs()[static_cast<std::size_t>(bit / IntegerBackend::limb_bits)] >>
             static_cast<unsigned>(bit % IntegerBackend::limb_bits)) &
            0xffU);
        hash.update(&byte, 1);
    }
}
void append_rational(boost::hash2::sha2_256& hash, const ExactRational& value) {
    const auto& b = detail::ExactFactory::value(value).backend();
    const unsigned char sign = b.num().sign() ? 1 : 0;
    hash.update(&sign, 1);
    append_integer(hash, b.num());
    append_integer(hash, b.denom());
}
PgGraphContentIdentity content_hash(const PgGraphData& d) {
    boost::hash2::sha2_256 hash;
    append_text(hash, "uvlas.pg.graph.v1");
    append_u64(hash, 1);
    append_rational(hash, d.x0);
    append_rational(hash, d.width);
    append_rational(hash, d.skew);
    append_u64(hash, d.height.terms().size());
    for (const auto& term : d.height.terms()) {
        append_u64(hash, static_cast<std::uint64_t>(term.powers.first));
        append_u64(hash, static_cast<std::uint64_t>(term.powers.second));
        append_rational(hash, term.coefficient);
    }
    append_u64(hash, static_cast<std::uint64_t>(d.n));
    append_u64(hash, d.cells.size());
    for (const auto& cell : d.cells.values()) {
        const auto a = static_cast<std::uint64_t>(cell.first),
                   b = static_cast<std::uint64_t>(cell.second),
                   axis = static_cast<std::uint64_t>(2 * d.n);
        for (const auto value : {a, b, a * axis + b, a, a + 1, b, b + 1, axis, 2 * axis})
            append_u64(hash, value);
    }
    append_text(hash, "anchor-s0-first-cell-v1");
    append_u64(hash, 0);
    const unsigned char negative = 1, one = 1;
    hash.update(&negative, 1);
    append_u64(hash, 1);
    hash.update(&one, 1);
    append_u64(hash, 1);
    hash.update(&one, 1);
    append_text(hash, "UnitSquare-v1");
    for (const auto value : {0U, 1U, 0U, 1U})
        append_u64(hash, value);
    const auto digest = hash.result();
    std::array<std::byte, 32> bytes{};
    for (std::size_t i = 0; i < bytes.size(); ++i)
        bytes[i] = static_cast<std::byte>(digest[i]);
    return PgGraphContentIdentity{bytes};
}
ExactRational small(std::int64_t numerator, std::uint64_t denominator, const PgGraphLimits& limits,
                    core::ResourcePlan& plan) {
    std::array<char, 32> n{}, d{};
    const auto nr = std::to_chars(n.data(), n.data() + n.size(), numerator),
               dr = std::to_chars(d.data(), d.data() + d.size(), denominator);
    if (nr.ec != std::errc{} || dr.ec != std::errc{})
        fail(core::FailureCode::ResourceExceeded);
    return take(ExactRational::from_integer_pair(
        {n.data(), static_cast<std::size_t>(nr.ptr - n.data())},
        {d.data(), static_cast<std::size_t>(dr.ptr - d.data())}, inner(limits), plan));
}
int compare_scaled(const IntegerBackend& a, std::uint64_t a_factor, const IntegerBackend& b,
                   std::uint64_t b_factor) noexcept {
    using Wide = boost::multiprecision::double_limb_type;
    using Limb = boost::multiprecision::limb_type;
    Wide acarry = 0, bcarry = 0;
    int comparison = 0;
    const auto size = std::max(a.size(), b.size());
    // Emit low-to-high with carry; the highest differing limb determines order.
    // No arbitrary-precision product or unadmitted scratch is constructed.
    for (std::size_t i = 0; i <= size; ++i) {
        const Wide av = (i < a.size() ? static_cast<Wide>(a.limbs()[i]) * a_factor : 0) + acarry;
        const Wide bv = (i < b.size() ? static_cast<Wide>(b.limbs()[i]) * b_factor : 0) + bcarry;
        const auto alow = static_cast<Limb>(av), blow = static_cast<Limb>(bv);
        if (alow != blow)
            comparison = alow < blow ? -1 : 1;
        acarry = av >> IntegerBackend::limb_bits;
        bcarry = bv >> IntegerBackend::limb_bits;
    }
    return comparison;
}
int compare_native(const ExactRational& value, core::Index numerator, core::Index denominator) {
    const auto sign = value.sign(), target = numerator < 0 ? -1 : numerator == 0 ? 0 : 1;
    if (sign != target)
        return sign < target ? -1 : 1;
    if (sign == 0)
        return 0;
    const auto& b = detail::ExactFactory::value(value).backend();
    const auto magnitude = numerator < 0 ? static_cast<std::uint64_t>(-(numerator + 1)) + 1
                                         : static_cast<std::uint64_t>(numerator);
    const auto result =
        compare_scaled(b.num(), static_cast<std::uint64_t>(denominator), b.denom(), magnitude);
    return sign < 0 ? -result : result;
}
core::Index locate_axis(const ExactRational& value, PgSide side, core::Index n) {
    check_side(side);
    if (compare_native(value, 0, 1) < 0 || compare_native(value, 1, 1) > 0)
        fail(core::FailureCode::InvalidInput);
    for (core::Index boundary = 0; boundary <= 2 * n; ++boundary) {
        const auto sign = compare_native(value, boundary, 2 * n);
        if (sign == 0) {
            if (side == PgSide::Lower && boundary > 0)
                return boundary - 1;
            if (side == PgSide::Upper && boundary < 2 * n)
                return boundary;
            fail(core::FailureCode::InvalidInput);
        }
        if (sign < 0) {
            if (side != PgSide::Interior)
                fail(core::FailureCode::InvalidInput);
            return boundary - 1;
        }
    }
    fail(core::FailureCode::InvalidInput);
}
void local_side(const ExactRational& value, PgSide side) {
    check_side(side);
    const auto lower = compare_native(value, -1, 1), upper = compare_native(value, 1, 1);
    if (lower < 0 || upper > 0)
        fail(core::FailureCode::InvalidInput);
    const auto expected = lower == 0   ? PgSide::Upper
                          : upper == 0 ? PgSide::Lower
                                       : PgSide::Interior;
    if (side != expected)
        fail(core::FailureCode::InvalidInput);
}
} // namespace

namespace detail {
const PgGraphData* PgGraphFactory::data(const PgGraph& owner) noexcept {
    return owner.data_.size() == 1 ? owner.data_.values().data() : nullptr;
}
PgGraph PgGraphFactory::graph(core::PlannedBuffer<PgGraphData> data,
                              PgGraphCallInfo info) noexcept {
    return PgGraph{std::move(data), info};
}
void PgGraphFactory::info(PgGraph& owner, PgGraphCallInfo info) noexcept {
    owner.info_ = info;
}
const PgCellPointData* PgGraphFactory::data(const PgCellPoint& owner) noexcept {
    return owner.data_.size() == 1 ? owner.data_.values().data() : nullptr;
}
PgCellPoint PgGraphFactory::point(core::PlannedBuffer<PgCellPointData> data,
                                  PgGraphCallInfo info) noexcept {
    return PgCellPoint{std::move(data), info};
}
void PgGraphFactory::info(PgCellPoint& owner, PgGraphCallInfo info) noexcept {
    owner.info_ = info;
}
const PgCellPolynomialsData* PgGraphFactory::data(const PgCellPolynomials& owner) noexcept {
    return owner.data_.size() == 1 ? owner.data_.values().data() : nullptr;
}
PgCellPolynomials PgGraphFactory::polynomials(core::PlannedBuffer<PgCellPolynomialsData> data,
                                              PgGraphCallInfo info) noexcept {
    return PgCellPolynomials{std::move(data), info};
}
void PgGraphFactory::info(PgCellPolynomials& owner, PgGraphCallInfo info) noexcept {
    owner.info_ = info;
}
const PgGraphValueData* PgGraphFactory::data(const PgGraphValue& owner) noexcept {
    return owner.data_.size() == 1 ? owner.data_.values().data() : nullptr;
}
PgGraphValue PgGraphFactory::value(core::PlannedBuffer<PgGraphValueData> data,
                                   PgGraphCallInfo info) noexcept {
    return PgGraphValue{std::move(data), info};
}
void PgGraphFactory::info(PgGraphValue& owner, PgGraphCallInfo info) noexcept {
    owner.info_ = info;
}
const PgRectangleActivityData*
PgGraphFactory::data(const QualifiedRectangularActivityMap& owner) noexcept {
    return owner.data_.size() == 1 ? owner.data_.values().data() : nullptr;
}
QualifiedRectangularActivityMap
PgGraphFactory::map(core::PlannedBuffer<PgRectangleActivityData> data,
                    PgGraphCallInfo info) noexcept {
    return QualifiedRectangularActivityMap{std::move(data), info};
}
void PgGraphFactory::info(QualifiedRectangularActivityMap& owner, PgGraphCallInfo info) noexcept {
    owner.info_ = info;
}
} // namespace detail
PgGraph::PgGraph(core::PlannedBuffer<detail::PgGraphData> data, PgGraphCallInfo info) noexcept
    : data_(std::move(data)), info_(info) {}
PgGraph::PgGraph(PgGraph&&) noexcept = default;
PgGraph& PgGraph::operator=(PgGraph&&) noexcept = default;
PgGraph::~PgGraph() = default;
const PgGraphCallInfo& PgGraph::preparation() const& noexcept {
    return info_;
}
PgCellPoint::PgCellPoint(core::PlannedBuffer<detail::PgCellPointData> data,
                         PgGraphCallInfo info) noexcept
    : data_(std::move(data)), info_(info) {}
PgCellPoint::PgCellPoint(PgCellPoint&&) noexcept = default;
PgCellPoint& PgCellPoint::operator=(PgCellPoint&&) noexcept = default;
PgCellPoint::~PgCellPoint() = default;
const PgGraphCallInfo& PgCellPoint::preparation() const& noexcept {
    return info_;
}
PgCellPolynomials::PgCellPolynomials(core::PlannedBuffer<detail::PgCellPolynomialsData> data,
                                     PgGraphCallInfo info) noexcept
    : data_(std::move(data)), info_(info) {}
PgCellPolynomials::PgCellPolynomials(PgCellPolynomials&&) noexcept = default;
PgCellPolynomials& PgCellPolynomials::operator=(PgCellPolynomials&&) noexcept = default;
PgCellPolynomials::~PgCellPolynomials() = default;
const PgGraphCallInfo& PgCellPolynomials::preparation() const& noexcept {
    return info_;
}
PgGraphValue::PgGraphValue(core::PlannedBuffer<detail::PgGraphValueData> data,
                           PgGraphCallInfo info) noexcept
    : data_(std::move(data)), info_(info) {}
PgGraphValue::PgGraphValue(PgGraphValue&&) noexcept = default;
PgGraphValue& PgGraphValue::operator=(PgGraphValue&&) noexcept = default;
PgGraphValue::~PgGraphValue() = default;
const PgGraphCallInfo& PgGraphValue::preparation() const& noexcept {
    return info_;
}
QualifiedRectangularActivityMap::QualifiedRectangularActivityMap(
    core::PlannedBuffer<detail::PgRectangleActivityData> data, PgGraphCallInfo info) noexcept
    : data_(std::move(data)), info_(info) {}
QualifiedRectangularActivityMap::QualifiedRectangularActivityMap(
    QualifiedRectangularActivityMap&&) noexcept = default;
QualifiedRectangularActivityMap&
QualifiedRectangularActivityMap::operator=(QualifiedRectangularActivityMap&&) noexcept = default;
QualifiedRectangularActivityMap::~QualifiedRectangularActivityMap() = default;
const PgGraphCallInfo& QualifiedRectangularActivityMap::preparation() const& noexcept {
    return info_;
}

core::Result<PgGraphView> PgGraph::view() const& noexcept {
    const auto* data = PgGraphFactory::data(*this);
    if (data == nullptr)
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, operation});
    const auto& id = data->identity;
    if (!id)
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, operation});
    return PgGraphView{data->x0, data->width,      data->skew,           data->height,
                       data->n,  data->convention, data->cells.values(), *id};
}
core::Result<PgCellPointView> PgCellPoint::view() const& noexcept {
    const auto* data = PgGraphFactory::data(*this);
    if (data == nullptr)
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, operation});
    const auto& snapshot = data->source;
    if (!snapshot)
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, operation});
    const auto definition = snapshot->view();
    if (!definition)
        return std::unexpected(definition.error());
    return PgCellPointView{data->s,    data->t,      data->u,      data->v,
                           data->cell, data->side_s, data->side_t, definition->identity};
}
core::Result<PgCellPolynomialsView> PgCellPolynomials::view() const& noexcept {
    const auto* data = PgGraphFactory::data(*this);
    if (data == nullptr)
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, operation});
    const auto& id = data->identity;
    if (!id)
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, operation});
    return PgCellPolynomialsView{data->components, data->cell, *id};
}
core::Result<PgGraphValueView> PgGraphValue::view() const& noexcept {
    const auto* data = PgGraphFactory::data(*this);
    if (data == nullptr)
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, operation});
    const auto& id = data->identity;
    if (!id)
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, operation});
    return PgGraphValueView{data->components, data->cell, data->orders, data->frame, *id};
}
core::Result<PgGraphView> QualifiedRectangularActivityMap::graph() const& noexcept {
    const auto* data = PgGraphFactory::data(*this);
    if (data == nullptr)
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, operation});
    const auto& snapshot = data->source;
    if (!snapshot)
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, operation});
    return snapshot->view();
}
core::Result<void>
QualifiedRectangularActivityMap::validate_for(const PgGraph& other) const noexcept {
    const auto* data = PgGraphFactory::data(*this);
    const auto* target = PgGraphFactory::data(other);
    if (data == nullptr || target == nullptr)
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, operation});
    const auto& snapshot = data->source;
    if (!snapshot)
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, operation});
    const auto* original = PgGraphFactory::data(*snapshot);
    if (original == nullptr)
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, operation});
    const auto& a = original->identity;
    const auto& b = target->identity;
    if (!a || !b)
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, operation});
    if (!(*a == *b) || !contents_equal(*original, *target))
        return std::unexpected(core::Failure{core::FailureCode::StaleGeneration, operation});
    return {};
}
core::Result<PgGraph> PgGraph::prepare(const PgGraphRequest& request, const PgGraphLimits& limits,
                                       core::ResourcePlan& plan) noexcept {
    return entry(limits, [&]() -> core::Result<PgGraph> {
        definition(request, limits);
        cancelled(limits);
        const auto before = plan.allocated_bytes();
        auto info = admit(prepare_bounds(request, limits), limits, plan);
        auto storage = take(plan.allocate<PgGraphData>(1, operation));
        auto& d = storage.values()[0];
        d.x0 = take(request.x0.clone(inner(limits), plan));
        cancelled(limits);
        d.width = take(request.width.clone(inner(limits), plan));
        d.skew = take(request.skew.clone(inner(limits), plan));
        d.height = take(request.height.clone(inner(limits), plan));
        cancelled(limits);
        d.n = request.n;
        d.convention = request.convention;
        const auto axis = 2 * d.n;
        d.cells = take(plan.allocate<PgCellIndex>(axis * axis, operation));
        for (core::Index a = 0; a < axis; ++a)
            for (core::Index b = 0; b < axis; ++b)
                d.cells.values()[static_cast<std::size_t>(a * axis + b)] = {a, b};
        d.identity = content_hash(d);
        cancelled(limits);
        return PgGraphFactory::graph(std::move(storage), completed(info, before, plan));
    });
}
core::Result<PgGraph> PgGraph::clone(const PgGraphLimits& limits,
                                     core::ResourcePlan& plan) const noexcept {
    const auto definition_view = view();
    if (!definition_view)
        return std::unexpected(definition_view.error());
    return prepare({definition_view->x0, definition_view->width, definition_view->skew,
                    definition_view->height, definition_view->n, definition_view->convention},
                   limits, plan);
}
core::Result<PgCellPoint> PgGraph::locate(const ExactRational& s, const ExactRational& t,
                                          PgSide side_s, PgSide side_t, const PgGraphLimits& limits,
                                          core::ResourcePlan& plan) const noexcept {
    return entry(limits, [&]() -> core::Result<PgCellPoint> {
        const auto& d = graph_data(*this, limits);
        check_value(s, limits);
        check_value(t, limits);
        check_side(side_s);
        check_side(side_t);
        cancelled(limits);
        const PgCellIndex cell{locate_axis(s, side_s, d.n), locate_axis(t, side_t, d.n)};
        const auto before = plan.allocated_bytes();
        auto info = admit(point_bounds(d, s, t, limits), limits, plan);
        auto storage = take(plan.allocate<PgCellPointData>(1, operation));
        auto& p = storage.values()[0];
        p.source = take(clone(nested(limits), plan));
        p.s = take(s.clone(inner(limits), plan));
        p.t = take(t.clone(inner(limits), plan));
        const auto scale = small(4 * d.n, 1, limits, plan),
                   cs = small(2 * cell.first + 1, 1, limits, plan),
                   ct = small(2 * cell.second + 1, 1, limits, plan);
        const auto ss = take(s.multiply(scale, inner(limits), plan)),
                   tt = take(t.multiply(scale, inner(limits), plan));
        p.u = take(ss.subtract(cs, inner(limits), plan));
        p.v = take(tt.subtract(ct, inner(limits), plan));
        p.cell = cell;
        p.side_s = side_s;
        p.side_t = side_t;
        cancelled(limits);
        return PgGraphFactory::point(std::move(storage), completed(info, before, plan));
    });
}
core::Result<PgCellPoint> PgGraph::from_local(PgCellIndex cell, const ExactRational& u,
                                              const ExactRational& v, PgSide side_s, PgSide side_t,
                                              const PgGraphLimits& limits,
                                              core::ResourcePlan& plan) const noexcept {
    return entry(limits, [&]() -> core::Result<PgCellPoint> {
        const auto& d = graph_data(*this, limits);
        check_cell(d, cell);
        check_value(u, limits);
        check_value(v, limits);
        local_side(u, side_s);
        local_side(v, side_t);
        cancelled(limits);
        const auto before = plan.allocated_bytes();
        auto info = admit(point_bounds(d, u, v, limits), limits, plan);
        auto storage = take(plan.allocate<PgCellPointData>(1, operation));
        auto& p = storage.values()[0];
        p.source = take(clone(nested(limits), plan));
        p.u = take(u.clone(inner(limits), plan));
        p.v = take(v.clone(inner(limits), plan));
        const auto scale = small(4 * d.n, 1, limits, plan),
                   cs = small(2 * cell.first + 1, 1, limits, plan),
                   ct = small(2 * cell.second + 1, 1, limits, plan);
        const auto us = take(u.add(cs, inner(limits), plan)),
                   vs = take(v.add(ct, inner(limits), plan));
        p.s = take(us.divide(scale, inner(limits), plan));
        p.t = take(vs.divide(scale, inner(limits), plan));
        p.cell = cell;
        p.side_s = side_s;
        p.side_t = side_t;
        cancelled(limits);
        return PgGraphFactory::point(std::move(storage), completed(info, before, plan));
    });
}
core::Result<PgCellPolynomials> PgGraph::cell_polynomials(PgCellIndex cell,
                                                          const PgGraphLimits& limits,
                                                          core::ResourcePlan& plan) const noexcept {
    return entry(limits, [&]() -> core::Result<PgCellPolynomials> {
        const auto& d = graph_data(*this, limits);
        check_cell(d, cell);
        cancelled(limits);
        const auto before = plan.allocated_bytes();
        auto info = admit(cell_bounds(d, limits), limits, plan);
        auto storage = take(plan.allocate<PgCellPolynomialsData>(1, operation));
        auto& output = storage.values()[0];
        const auto scale = small(1, static_cast<std::uint64_t>(4 * d.n), limits, plan);
        const auto cs = small(2 * cell.first + 1, 1, limits, plan),
                   ct = small(2 * cell.second + 1, 1, limits, plan);
        const auto os = take(cs.multiply(scale, inner(limits), plan)),
                   ot = take(ct.multiply(scale, inner(limits), plan));
        const std::array st{ExactInputTerm{{0, 0}, os}, ExactInputTerm{{1, 0}, scale}};
        const std::array tt{ExactInputTerm{{0, 0}, ot}, ExactInputTerm{{0, 1}, scale}};
        const auto s = take(ExactCellPolynomial::make(st, inner(limits), plan)),
                   t = take(ExactCellPolynomial::make(tt, inner(limits), plan));
        const auto su = take(d.width.multiply(scale, inner(limits), plan)),
                   tv = take(d.skew.multiply(scale, inner(limits), plan));
        const auto wo = take(d.width.multiply(os, inner(limits), plan)),
                   ao = take(d.skew.multiply(ot, inner(limits), plan));
        const auto wa = take(wo.add(ao, inner(limits), plan)),
                   constant = take(wa.add(d.x0, inner(limits), plan));
        const std::array xt{ExactInputTerm{{0, 0}, constant}, ExactInputTerm{{1, 0}, su},
                            ExactInputTerm{{0, 1}, tv}};
        output.components[0] = take(ExactCellPolynomial::make(xt, inner(limits), plan));
        output.components[1] = take(t.clone(inner(limits), plan));
        cancelled(limits);
        output.components[2] = take(d.height.substitute(s, t, inner(limits), plan));
        output.cell = cell;
        output.identity = identity(d);
        cancelled(limits);
        return PgGraphFactory::polynomials(std::move(storage), completed(info, before, plan));
    });
}
core::Result<PgGraphValue> PgGraph::evaluate(const PgCellPoint& point, PgDerivativeOrders orders,
                                             PgDerivativeFrame frame, const PgGraphLimits& limits,
                                             core::ResourcePlan& plan) const noexcept {
    return entry(limits, [&]() -> core::Result<PgGraphValue> {
        const auto& d = graph_data(*this, limits);
        check_orders(orders, frame);
        const auto* p = PgGraphFactory::data(point);
        if (p == nullptr)
            fail(core::FailureCode::InvalidInput);
        const auto& snapshot = p->source;
        if (!snapshot)
            fail(core::FailureCode::InvalidInput);
        const auto& source = graph_data(*snapshot, limits);
        if (!(identity(d) == identity(source)) || !contents_equal(d, source))
            fail(core::FailureCode::StaleGeneration);
        for (const auto* coordinate : {&p->s, &p->t, &p->u, &p->v})
            check_value(*coordinate, limits);
        cancelled(limits);
        const auto before = plan.allocated_bytes();
        auto info = admit(evaluation_bounds(d, *p, orders, frame, limits), limits, plan);
        auto storage = take(plan.allocate<PgGraphValueData>(1, operation));
        auto& value = storage.values()[0];
        value.cell = p->cell;
        value.orders = orders;
        value.frame = frame;
        value.identity = identity(d);
        if (!structural_zero(d, orders)) {
            const auto total = orders.first + orders.second;
            if (total == 0) {
                const auto ws = take(d.width.multiply(p->s, inner(limits), plan)),
                           at = take(d.skew.multiply(p->t, inner(limits), plan));
                const auto sum = take(ws.add(at, inner(limits), plan));
                value.components[0] = take(sum.add(d.x0, inner(limits), plan));
                value.components[1] = take(p->t.clone(inner(limits), plan));
            } else if (orders.first == 1 && orders.second == 0)
                value.components[0] = take(d.width.clone(inner(limits), plan));
            else if (orders.first == 0 && orders.second == 1) {
                value.components[0] = take(d.skew.clone(inner(limits), plan));
                value.components[1] = small(1, 1, limits, plan);
            }
            const ExactChartPolynomial* derivative = &d.height;
            std::optional<ExactChartPolynomial> first, second;
            if (orders.first != 0) {
                first = take(
                    derivative->differentiate(ExactAxis::First, orders.first, inner(limits), plan));
                derivative = &*first;
            }
            if (orders.second != 0) {
                second = take(derivative->differentiate(ExactAxis::Second, orders.second,
                                                        inner(limits), plan));
                derivative = &*second;
            }
            if (!derivative->terms().empty())
                value.components[2] = take(derivative->evaluate(p->s, p->t, inner(limits), plan));
            const auto nonzero =
                std::any_of(value.components.begin(), value.components.end(),
                            [](const ExactRational& item) { return item.sign() != 0; });
            // A nonzero derivative has total<=4. Never form this power for a zero jet.
            if (frame == PgDerivativeFrame::CellUV && total != 0 && nonzero) {
                std::uint64_t denominator = 1;
                for (core::Index k = 0; k < total; ++k)
                    denominator = mul(denominator, static_cast<std::uint64_t>(4 * d.n));
                const auto scale = small(1, denominator, limits, plan);
                for (auto& component : value.components)
                    if (component.sign() != 0)
                        component = take(component.multiply(scale, inner(limits), plan));
            }
        }
        cancelled(limits);
        return PgGraphFactory::value(std::move(storage), completed(info, before, plan));
    });
}
core::Result<PgGraphValue> PgGraph::anchor_at(const ExactRational& t, PgSide side_t,
                                              core::Index order, const PgGraphLimits& limits,
                                              core::ResourcePlan& plan) const noexcept {
    return entry(limits, [&]() -> core::Result<PgGraphValue> {
        const auto& d = graph_data(*this, limits);
        check_orders({0, order}, PgDerivativeFrame::ChartST);
        check_value(t, limits);
        (void)locate_axis(t, side_t, d.n);
        cancelled(limits);
        const ExactRational zero;
        const auto pb = point_bounds(d, zero, t, limits);
        const auto m = d.height.terms().size();
        const auto f = maximum_bits(d.height, limits), v = graph_bits(d, limits),
                   p = add(value_bits(t), 7);
        const auto ew =
            structural_zero(d, {0, order})
                ? 0
                : std::max({mul(m, add(f, 13)), add(mul(m, add(add(f, 12), add(mul(4, p), 1))), 17),
                            add(add(mul(3, v), mul(2, p)), 20), std::uint64_t{21}});
        const auto ep =
            structural_zero(d, {0, order})
                ? sizeof(PgGraphValueData)
                : add(add(sizeof(PgGraphValueData), mul(8, sizeof(detail::RationalStorage))),
                      mul(2, poly_payload(m)));
        const Bounds bounds{std::max(pb.work, ew), std::max(pb.raw, m), 64, add(pb.payload, ep),
                            add(retained(d), retained(t))};
        const auto before = plan.allocated_bytes();
        auto info = admit(bounds, limits, plan);
        const auto point = take(locate(zero, t, PgSide::Upper, side_t, nested(limits), plan));
        auto result =
            take(evaluate(point, {0, order}, PgDerivativeFrame::ChartST, nested(limits), plan));
        cancelled(limits);
        PgGraphFactory::info(result, completed(info, before, plan));
        return result;
    });
}
core::Result<QualifiedRectangularActivityMap>
QualifiedRectangularActivityMap::prepare(const PgGraph& graph, const PgGraphLimits& limits,
                                         core::ResourcePlan& plan) noexcept {
    return entry(limits, [&]() -> core::Result<QualifiedRectangularActivityMap> {
        const auto& d = graph_data(graph, limits);
        cancelled(limits);
        native_chars(d.cells.size(), limits);
        const auto base =
            prepare_bounds({d.x0, d.width, d.skew, d.height, d.n, d.convention}, limits);
        // Native count<=64 uses at most two decimal digits (4*2+1=9); area*count uses <=8 bits.
        const Bounds bounds{std::max(base.work, std::uint64_t{9}), base.raw, 22,
                            add(add(prepare_payload(d), sizeof(PgRectangleActivityData)),
                                mul(4, sizeof(detail::RationalStorage))),
                            retained(d)};
        const auto before = plan.allocated_bytes();
        auto info = admit(bounds, limits, plan);
        auto storage = take(plan.allocate<PgRectangleActivityData>(1, operation));
        auto& output = storage.values()[0];
        output.source = take(graph.clone(nested(limits), plan));
        const auto count = static_cast<std::uint64_t>(d.cells.size());
        output.area = small(1, count, limits, plan);
        const auto count_value = small(static_cast<std::int64_t>(count), 1, limits, plan),
                   one = small(1, 1, limits, plan);
        output.total = take(output.area.multiply(count_value, inner(limits), plan));
        if (output.area.sign() <= 0 || !(output.total == one))
            fail(core::FailureCode::InvalidInput);
        cancelled(limits);
        return PgGraphFactory::map(std::move(storage), completed(info, before, plan));
    });
}
namespace {
core::Result<ExactRational> area_copy(const QualifiedRectangularActivityMap& map, bool total,
                                      const PgGraphLimits& limits,
                                      core::ResourcePlan& plan) noexcept {
    return entry(limits, [&]() -> core::Result<ExactRational> {
        const auto* data = PgGraphFactory::data(map);
        if (data == nullptr)
            fail(core::FailureCode::InvalidInput);
        const auto& snapshot = data->source;
        if (!snapshot)
            fail(core::FailureCode::InvalidInput);
        const auto& graph = graph_data(*snapshot, limits);
        check_value(data->area, limits);
        check_value(data->total, limits);
        const auto& value = total ? data->total : data->area;
        cancelled(limits);
        const Bounds bounds{value_bits(value),
                            0,
                            1,
                            sizeof(detail::RationalStorage),
                            add(retained(graph), add(retained(data->area), retained(data->total))),
                            true};
        const auto info = admit(bounds, limits, plan);
        auto result =
            detail::ExactFactory::rational(detail::Rational{detail::ExactFactory::value(value)},
                                           inner(limits), plan, info.arithmetic);
        cancelled(limits);
        return result;
    });
}
} // namespace
core::Result<ExactRational>
QualifiedRectangularActivityMap::cell_area(const PgGraphLimits& limits,
                                           core::ResourcePlan& plan) const noexcept {
    return area_copy(*this, false, limits, plan);
}
core::Result<ExactRational>
QualifiedRectangularActivityMap::total_area(const PgGraphLimits& limits,
                                            core::ResourcePlan& plan) const noexcept {
    return area_copy(*this, true, limits, plan);
}
} // namespace uvlas::structure
