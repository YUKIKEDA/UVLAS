#include "exact_polynomial_detail.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <limits>
#include <type_traits>
#include <uvlas/structure/exact_polynomial.hpp>

namespace uvlas::structure {
namespace {
constexpr auto operation = core::Operation::PrepareModel;
using detail::ExactFactory;
using detail::Integer;
using detail::IntegerBackend;
using detail::Rational;
struct Stop {
    core::Failure failure;
};
[[noreturn]] void stop(core::FailureCode code) {
    throw Stop{core::Failure{code, operation}};
}
std::uint64_t plus(std::uint64_t a, std::uint64_t b) {
    if (b > std::numeric_limits<std::uint64_t>::max() - a)
        stop(core::FailureCode::ResourceExceeded);
    return a + b;
}
std::uint64_t times(std::uint64_t a, std::uint64_t b) {
    if (a != 0 && b > std::numeric_limits<std::uint64_t>::max() / a)
        stop(core::FailureCode::ResourceExceeded);
    return a * b;
}
core::Index index(std::uint64_t value) {
    if (value > static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max()))
        stop(core::FailureCode::ResourceExceeded);
    return static_cast<core::Index>(value);
}
template <class T>
T take(core::Result<T> value) {
    if (!value)
        throw Stop{value.error()};
    return std::move(*value);
}
std::uint64_t bits(const IntegerBackend& value) noexcept {
    if (value.size() == 1 && value.limbs()[0] == 0)
        return 0;
    return static_cast<std::uint64_t>(
        (value.size() - 1) * IntegerBackend::limb_bits + IntegerBackend::limb_bits -
        static_cast<std::size_t>(std::countl_zero(value.limbs()[value.size() - 1])));
}
std::uint64_t value_bits(const Rational& value) noexcept {
    return std::max(bits(value.backend().num()), bits(value.backend().denom()));
}
bool zero(const Rational& value) noexcept {
    return bits(value.backend().num()) == 0;
}
void canonical_zero(Rational& value) {
    // Pinned Boost mixed-backend division can retain 0/k. Never publish it.
    if (zero(value))
        value.backend().denom() = Integer{1}.backend();
}
std::uint64_t retained_bytes(const Rational& value) {
    const auto capacity = [](const IntegerBackend& v) {
        return v.capacity() > IntegerBackend::internal_limb_count ? v.capacity() : 0;
    };
    return times(plus(capacity(value.backend().num()), capacity(value.backend().denom())),
                 sizeof(boost::multiprecision::limb_type));
}
std::uint64_t retained_bytes(std::span<const ExactTerm> terms) {
    std::uint64_t result = 0;
    for (const auto& term : terms)
        result = plus(result, retained_bytes(ExactFactory::value(term.coefficient)));
    return result;
}
void validate(const ExactLimits& limits) {
    if (limits.external_memory.valueless_by_exception() ||
        (limits.timing != ExactTimingMode::Disabled && limits.timing != ExactTimingMode::Enabled))
        stop(core::FailureCode::InvalidInput);
    (void)index(limits.maximum_terms);
    (void)index(limits.maximum_raw_terms);
    (void)index(limits.maximum_axis_degree);
    (void)index(limits.maximum_total_degree);
    (void)index(limits.maximum_input_characters);
    (void)index(limits.maximum_value_bits);
    (void)index(limits.maximum_work_bits);
}
void check_value(const Rational& value, const ExactLimits& limits) {
    if (value_bits(value) > limits.maximum_value_bits)
        stop(core::FailureCode::ResourceExceeded);
}
void check_powers(ExactPowers powers, const ExactLimits& limits) {
    if (powers.first < 0 || powers.second < 0)
        stop(core::FailureCode::InvalidInput);
    const auto first = static_cast<std::uint64_t>(powers.first),
               second = static_cast<std::uint64_t>(powers.second);
    if (first > limits.maximum_axis_degree || second > limits.maximum_axis_degree ||
        plus(first, second) > limits.maximum_total_degree)
        stop(core::FailureCode::ResourceExceeded);
}
void check_axis(ExactAxis axis) {
    if (axis != ExactAxis::First && axis != ExactAxis::Second)
        stop(core::FailureCode::InvalidInput);
}
enum class Scratch : std::uint8_t { Rational, Normalize, Compare, Copy };
ExactPreparationInfo admit(const ExactLimits& limits, std::uint64_t work, std::uint64_t raw,
                           std::uint64_t borrowed, std::uint64_t slots = 7,
                           Scratch scratch = Scratch::Rational) {
    if (work > limits.maximum_work_bits || raw > limits.maximum_raw_terms)
        stop(core::FailureCode::ResourceExceeded);
    if (std::holds_alternative<ExactRequireProvenBound>(limits.external_memory))
        stop(core::FailureCode::Unsupported);
    constexpr auto width = IntegerBackend::limb_bits;
    const auto limbs =
        slots == 0 && work == 0 ? 0 : plus(plus(work / width, work % width != 0 ? 1 : 0), 2);
    const auto guarded_bits = times(limbs, width);
    if (guarded_bits > limits.maximum_work_bits)
        stop(core::FailureCode::ResourceExceeded);
    const auto capacity = limbs <= IntegerBackend::internal_limb_count ? 0 : times(4, limbs);
    const auto old_heap = capacity == 0 ? 0 : limbs;
    // GCD's shared aliases can become owners after fallback swaps. Count U/V/t
    // as well as the fallback and long-division temporaries, not only the pool.
    const auto gcd = plus(times(6, plus(limbs, 1)), times(5, capacity));
    const auto division = times(2, capacity);
    const auto karatsuba = limbs >= 60 ? times(5, limbs) : 0;
    const auto nested = scratch == Scratch::Copy        ? 0
                        : scratch == Scratch::Compare   ? karatsuba
                        : scratch == Scratch::Normalize ? std::max(gcd, division)
                                                        : std::max({gcd, division, karatsuba});
    const auto total = plus(borrowed, times(plus(plus(times(slots, capacity), nested), old_heap),
                                            sizeof(boost::multiprecision::limb_type)));
    const auto* permission = std::get_if<ExactEstimateAdmission>(&limits.external_memory);
    if (permission == nullptr)
        stop(core::FailureCode::InvalidInput);
    if (total > permission->bytes)
        stop(core::FailureCode::ResourceExceeded);
    return {total, work, guarded_bits, raw, std::nullopt};
}
template <class F>
auto entry(const ExactLimits& limits, F&& action) noexcept -> std::invoke_result_t<F> {
    using Result = std::invoke_result_t<F>;
    core::StageRecorder recorder{limits.timing == ExactTimingMode::Enabled, operation};
    auto result =
        core::guard_entry(operation, core::FailureCode::InvalidInput, {}, [&]() -> Result {
            try {
                validate(limits);
                const core::StageScope scope{recorder, core::MeasurementStage::Preparation};
                return action();
            } catch (const Stop& failure) {
                return std::unexpected(failure.failure);
            }
        });
    if (result) {
        const auto recording = recorder.snapshot();
        if (!recording)
            return std::unexpected(recording.error());
        if constexpr (requires {
                          ExactFactory::timing(*result, std::optional<core::StageTiming>{});
                      }) {
            ExactFactory::timing(*result, recorder.enabled()
                                              ? std::optional{(*recording)[static_cast<std::size_t>(
                                                    core::MeasurementStage::Preparation)]}
                                              : std::nullopt);
        }
    }
    return result;
}
enum class Arithmetic : std::uint8_t { Add, Subtract, Multiply, Divide };
std::uint64_t rational_bound(const Rational& a, const Rational& b, Arithmetic kind) {
    const auto an = bits(a.backend().num()), ad = bits(a.backend().denom());
    const auto bn = bits(b.backend().num()), bd = bits(b.backend().denom());
    std::uint64_t n = 0, d = 0;
    if (kind == Arithmetic::Add || kind == Arithmetic::Subtract) {
        n = plus(std::max(plus(an, bd), plus(bn, ad)), 1);
        d = plus(ad, bd);
    } else if (kind == Arithmetic::Multiply) {
        n = plus(an, bn);
        d = plus(ad, bd);
    } else {
        n = plus(an, bd);
        d = plus(ad, bn);
    }
    return std::max({n, d, an, ad, bn, bd});
}
Rational calculate(const Rational& a, const Rational& b, Arithmetic kind) {
    Rational result;
    // Fresh nonalias results and et_off are part of the reviewed live-slot model.
    switch (kind) {
    case Arithmetic::Add:
        boost::multiprecision::backends::eval_add(result.backend(), a.backend(), b.backend());
        break;
    case Arithmetic::Subtract:
        boost::multiprecision::backends::eval_subtract(result.backend(), a.backend(), b.backend());
        break;
    case Arithmetic::Multiply:
        boost::multiprecision::backends::eval_multiply(result.backend(), a.backend(), b.backend());
        break;
    case Arithmetic::Divide:
        if (zero(b))
            stop(core::FailureCode::InvalidInput);
        boost::multiprecision::backends::eval_divide(result.backend(), a.backend(), b.backend());
        break;
    }
    canonical_zero(result);
    return result;
}
struct Digits {
    std::string_view magnitude;
    bool negative;
};
Digits integer_digits(std::string_view text, const ExactLimits& limits) {
    if (text.size() > limits.maximum_input_characters)
        stop(core::FailureCode::ResourceExceeded);
    if (text.empty())
        stop(core::FailureCode::InvalidInput);
    bool negative = false;
    if (text.front() == '+' || text.front() == '-') {
        negative = text.front() == '-';
        text.remove_prefix(1);
    }
    if (text.empty())
        stop(core::FailureCode::InvalidInput);
    for (const char c : text)
        if (c < '0' || c > '9')
            stop(core::FailureCode::InvalidInput);
    while (text.size() > 1 && text.front() == '0')
        text.remove_prefix(1);
    if (text == "0")
        negative = false;
    return {text, negative};
}
Integer parse_integer(Digits input) {
    // Decimal block accumulation avoids cpp_int's leading-zero octal grammar and
    // needs no NUL-terminated copy of the borrowed input.
    Integer result = 0;
    constexpr std::size_t block_digits = boost::multiprecision::digits_per_block_10;
    for (std::size_t start = 0; start < input.magnitude.size();) {
        const auto count = std::min(block_digits, input.magnitude.size() - start);
        boost::multiprecision::limb_type block = 0, multiplier = 1;
        for (std::size_t i = 0; i < count; ++i) {
            block = block * 10 + static_cast<unsigned>(input.magnitude[start + i] - '0');
            multiplier *= 10;
        }
        result *= multiplier;
        result += block;
        start += count;
    }
    if (input.negative)
        result.backend().negate();
    return result;
}
Rational fraction(Integer numerator, Integer denominator) {
    if (bits(denominator.backend()) == 0)
        stop(core::FailureCode::InvalidInput);
    Rational result;
    result.backend().num() = std::move(numerator.backend());
    result.backend().denom() = std::move(denominator.backend());
    result.backend().normalize();
    canonical_zero(result);
    return result;
}
Integer integer_power(Integer base, std::uint64_t exponent) {
    Integer result = 1;
    while (exponent != 0) {
        if ((exponent & 1) != 0) {
            Integer product;
            boost::multiprecision::backends::eval_multiply(product.backend(), result.backend(),
                                                           base.backend());
            result = std::move(product);
        }
        exponent >>= 1;
        if (exponent != 0) {
            Integer square;
            boost::multiprecision::backends::eval_multiply(square.backend(), base.backend(),
                                                           base.backend());
            base = std::move(square);
        }
    }
    return result;
}
Rational rational_power(const Rational& base, std::uint64_t exponent) {
    Rational result = 1, power = base;
    while (exponent != 0) {
        if ((exponent & 1) != 0)
            result = calculate(result, power, Arithmetic::Multiply);
        exponent >>= 1;
        if (exponent != 0)
            power = calculate(power, power, Arithmetic::Multiply);
    }
    return result;
}
struct RawTerm {
    ExactPowers powers{0, 0};
    Rational value;
};
using Raw = core::PlannedBuffer<RawTerm>;
Raw raw(std::uint64_t size, core::ResourcePlan& plan) {
    return take(plan.allocate<RawTerm>(index(size), operation));
}
bool before(const RawTerm& a, const RawTerm& b) noexcept {
    return a.powers.first < b.powers.first ||
           (a.powers.first == b.powers.first && a.powers.second < b.powers.second);
}
std::uint64_t maximum_bits(std::span<const ExactTerm> terms, const ExactLimits& limits) {
    std::uint64_t result = 1;
    for (const auto& term : terms) {
        check_powers(term.powers, limits);
        const auto& v = ExactFactory::value(term.coefficient);
        check_value(v, limits);
        result = std::max(result, value_bits(v));
    }
    return result;
}
std::uint64_t merge_bound(std::uint64_t raw_count, std::uint64_t coefficient_bound) {
    return raw_count == 0 ? 0 : times(raw_count, plus(coefficient_bound, 1));
}
core::Result<detail::PolynomialValue> finish(Raw values, const ExactLimits& limits,
                                             core::ResourcePlan& plan, ExactPreparationInfo info) {
    auto data = values.values();
    if (data.size() > 1)
        std::sort(data.begin(), data.end(), before);
    std::size_t used = 0;
    for (std::size_t i = 0; i < data.size(); ++i) {
        if (used != 0 && data[used - 1].powers.first == data[i].powers.first &&
            data[used - 1].powers.second == data[i].powers.second) {
            data[used - 1].value = calculate(data[used - 1].value, data[i].value, Arithmetic::Add);
        } else {
            if (used != i)
                data[used] = std::move(data[i]);
            ++used;
        }
    }
    std::size_t count = 0;
    for (std::size_t i = 0; i < used; ++i)
        if (!zero(data[i].value)) {
            check_value(data[i].value, limits);
            ++count;
        }
    if (count > limits.maximum_terms)
        stop(core::FailureCode::ResourceExceeded);
    auto output = take(plan.allocate<ExactTerm>(index(count), operation));
    std::size_t next = 0;
    for (std::size_t i = 0; i < used; ++i)
        if (!zero(data[i].value)) {
            output.values()[next].powers = data[i].powers;
            output.values()[next].coefficient =
                take(ExactFactory::rational(std::move(data[i].value), limits, plan, info));
            ++next;
        }
    return ExactFactory::polynomial(std::move(output), plan, info);
}
ExactPreparationInfo poly_plan(const ExactLimits& limits, std::uint64_t work, std::uint64_t count,
                               std::uint64_t borrowed) {
    // Substitution can retain the global raw output, the current expansion and
    // its replacement: three arrays, two integers per coefficient => 6*count.
    // Evaluation/powers retain accumulator, first/second powers, power result/base
    // (ten integers); one fresh rational multiply adds its two output components
    // and four helper backends (six), giving the fixed sixteen-slot allowance.
    // Sorting/merging/publication use no larger live set. Borrowed capacities are
    // counted separately and the nested GCD/divide/Karatsuba scratch is a max.
    if (count == 0)
        return admit(limits, 0, 0, borrowed, 0, Scratch::Copy);
    return admit(limits, work, count, borrowed, plus(16, times(6, count)));
}
} // namespace

namespace detail {
const Rational& ExactFactory::value(const ExactRational& value) {
    static const Rational zero_value;
    return value.storage_.size() == 0 ? zero_value : value.storage_.values()[0].value;
}
core::Result<ExactRational> ExactFactory::rational(Rational value, const ExactLimits& limits,
                                                   core::ResourcePlan& plan,
                                                   ExactPreparationInfo info) {
    check_value(value, limits);
    canonical_zero(value);
    if (zero(value))
        return ExactRational{core::PlannedBuffer<RationalStorage>{}, info};
    auto storage = take(plan.allocate<RationalStorage>(1, operation));
    storage.values()[0].value = std::move(value);
    return ExactRational{std::move(storage), info};
}
core::Result<PolynomialValue> ExactFactory::polynomial(core::PlannedBuffer<ExactTerm> terms,
                                                       core::ResourcePlan& plan,
                                                       ExactPreparationInfo info) {
    if (terms.size() == 0)
        return PolynomialValue{core::PlannedBuffer<PolynomialStorage>{}, info};
    auto storage = take(plan.allocate<PolynomialStorage>(1, operation));
    storage.values()[0].terms = std::move(terms);
    return PolynomialValue{std::move(storage), info};
}
void ExactFactory::timing(ExactRational& value, std::optional<core::StageTiming> timing) noexcept {
    value.info_.timing = timing;
}
void ExactFactory::timing(PolynomialValue& value,
                          std::optional<core::StageTiming> timing) noexcept {
    value.info_.timing = timing;
}
PolynomialValue::PolynomialValue() noexcept = default;
PolynomialValue::PolynomialValue(core::PlannedBuffer<PolynomialStorage> storage,
                                 ExactPreparationInfo info) noexcept
    : storage_(std::move(storage)), info_(info) {}
PolynomialValue::PolynomialValue(PolynomialValue&&) noexcept = default;
PolynomialValue& PolynomialValue::operator=(PolynomialValue&&) noexcept = default;
PolynomialValue::~PolynomialValue() = default;
std::span<const ExactTerm> PolynomialValue::terms() const& noexcept {
    return storage_.size() == 0 ? std::span<const ExactTerm>{}
                                : storage_.values()[0].terms.values();
}
const ExactPreparationInfo& PolynomialValue::preparation() const& noexcept {
    return info_;
}
} // namespace detail

ExactRational::ExactRational() noexcept = default;
ExactRational::ExactRational(core::PlannedBuffer<detail::RationalStorage> storage,
                             ExactPreparationInfo info) noexcept
    : storage_(std::move(storage)), info_(info) {}
ExactRational::ExactRational(ExactRational&&) noexcept = default;
ExactRational& ExactRational::operator=(ExactRational&&) noexcept = default;
ExactRational::~ExactRational() = default;
const ExactPreparationInfo& ExactRational::preparation() const& noexcept {
    return info_;
}
int ExactRational::sign() const noexcept {
    if (storage_.size() == 0)
        return 0;
    return storage_.values()[0].value.backend().num().sign() ? -1 : 1;
}
std::uint64_t ExactRational::numerator_bits() const noexcept {
    return storage_.size() == 0 ? 0 : bits(storage_.values()[0].value.backend().num());
}
std::uint64_t ExactRational::denominator_bits() const noexcept {
    return storage_.size() == 0 ? 1 : bits(storage_.values()[0].value.backend().denom());
}
bool ExactRational::operator==(const ExactRational& other) const noexcept {
    if (storage_.size() == 0 || other.storage_.size() == 0)
        return storage_.size() == other.storage_.size();
    const auto& a = storage_.values()[0].value;
    const auto& b = other.storage_.values()[0].value;
    return a.backend().num().compare(b.backend().num()) == 0 &&
           a.backend().denom().compare(b.backend().denom()) == 0;
}
core::Result<ExactRational> ExactRational::from_integer_pair(std::string_view numerator,
                                                             std::string_view denominator,
                                                             const ExactLimits& limits,
                                                             core::ResourcePlan& plan) noexcept {
    return entry(limits, [&]() -> core::Result<ExactRational> {
        const auto n = integer_digits(numerator, limits), d = integer_digits(denominator, limits);
        if (d.magnitude == "0")
            stop(core::FailureCode::InvalidInput);
        const auto work = plus(times(4, std::max(n.magnitude.size(), d.magnitude.size())), 1);
        const auto info = admit(limits, work, 0, 0, 6, Scratch::Normalize);
        return ExactFactory::rational(fraction(parse_integer(n), parse_integer(d)), limits, plan,
                                      info);
    });
}
core::Result<ExactRational> ExactRational::from_binary64(double input, const ExactLimits& limits,
                                                         core::ResourcePlan& plan) noexcept {
    return entry(limits, [&]() -> core::Result<ExactRational> {
        const auto word = std::bit_cast<std::uint64_t>(input);
        const auto exponent = (word >> 52) & 0x7ff;
        const auto fraction_bits = word & UINT64_C(0x000fffffffffffff);
        if (exponent == 0x7ff)
            stop(core::FailureCode::Nonfinite);
        const auto power = exponent == 0 ? -1074 : static_cast<int>(exponent) - 1075;
        const auto mantissa =
            exponent == 0 ? fraction_bits : fraction_bits | UINT64_C(0x0010000000000000);
        const auto work =
            mantissa == 0
                ? 1
                : std::max<std::uint64_t>(53 + (power > 0 ? static_cast<unsigned>(power) : 0),
                                          power < 0 ? static_cast<unsigned>(-power) + 1 : 1);
        const auto info = admit(limits, work, 0, 0, 6);
        Integer n = mantissa, d = 1;
        if (mantissa != 0) {
            if (power > 0)
                n <<= static_cast<unsigned>(power);
            else if (power < 0)
                d <<= static_cast<unsigned>(-power);
            if ((word >> 63) != 0)
                n.backend().negate();
        }
        return ExactFactory::rational(fraction(std::move(n), std::move(d)), limits, plan, info);
    });
}
core::Result<ExactRational> ExactRational::from_decimal(std::string_view input,
                                                        const ExactLimits& limits,
                                                        core::ResourcePlan& plan) noexcept {
    return entry(limits, [&]() -> core::Result<ExactRational> {
        if (input.size() > limits.maximum_input_characters)
            stop(core::FailureCode::ResourceExceeded);
        if (input.empty())
            stop(core::FailureCode::InvalidInput);
        bool negative = false;
        if (input.front() == '+' || input.front() == '-') {
            negative = input.front() == '-';
            input.remove_prefix(1);
        }
        std::size_t position = 0, digits = 0, fractional = 0;
        bool dot = false, nonzero = false;
        while (position < input.size() && input[position] != 'e' && input[position] != 'E') {
            const auto c = input[position++];
            if (c == '.' && !dot) {
                dot = true;
                continue;
            }
            if (c < '0' || c > '9')
                stop(core::FailureCode::InvalidInput);
            ++digits;
            if (dot)
                ++fractional;
            nonzero = nonzero || c != '0';
        }
        if (digits == 0)
            stop(core::FailureCode::InvalidInput);
        const auto mantissa_end = position;
        bool negative_exponent = false;
        std::uint64_t exponent = 0;
        if (position < input.size()) {
            ++position;
            if (position < input.size() && (input[position] == '+' || input[position] == '-')) {
                negative_exponent = input[position] == '-';
                ++position;
            }
            if (position == input.size())
                stop(core::FailureCode::InvalidInput);
            for (; position < input.size(); ++position) {
                const auto c = input[position];
                if (c < '0' || c > '9')
                    stop(core::FailureCode::InvalidInput);
                exponent = plus(times(exponent, 10), static_cast<unsigned>(c - '0'));
            }
        }
        const auto scale = negative_exponent        ? plus(exponent, fractional)
                           : exponent >= fractional ? exponent - fractional
                                                    : fractional - exponent;
        const bool divide_power = negative_exponent || exponent < fractional;
        const auto work = nonzero ? plus(times(4, plus(digits, scale)), 1) : 1;
        const auto info = admit(limits, work, 0, 0, 6);
        Integer coefficient = 0;
        for (std::size_t i = 0; i < mantissa_end; ++i)
            if (input[i] != '.') {
                coefficient *= 10;
                coefficient += static_cast<unsigned>(input[i] - '0');
            }
        Integer denominator = 1;
        if (nonzero && scale != 0) {
            auto power = integer_power(Integer{10}, scale);
            if (divide_power)
                denominator = std::move(power);
            else {
                Integer product;
                boost::multiprecision::backends::eval_multiply(
                    product.backend(), coefficient.backend(), power.backend());
                coefficient = std::move(product);
            }
        }
        if (negative)
            coefficient.backend().negate();
        return ExactFactory::rational(fraction(std::move(coefficient), std::move(denominator)),
                                      limits, plan, info);
    });
}
namespace {
core::Result<ExactRational> binary_rational(const ExactRational& left, const ExactRational& right,
                                            Arithmetic kind, const ExactLimits& limits,
                                            core::ResourcePlan& plan) noexcept {
    return entry(limits, [&]() -> core::Result<ExactRational> {
        const auto& a = ExactFactory::value(left);
        const auto& b = ExactFactory::value(right);
        check_value(a, limits);
        check_value(b, limits);
        if (kind == Arithmetic::Divide && zero(b))
            stop(core::FailureCode::InvalidInput);
        const auto info =
            admit(limits, rational_bound(a, b, kind), 0, plus(retained_bytes(a), retained_bytes(b)),
                  kind == Arithmetic::Add || kind == Arithmetic::Subtract ? 7 : 6);
        return ExactFactory::rational(calculate(a, b, kind), limits, plan, info);
    });
}
} // namespace
core::Result<ExactRational> ExactRational::add(const ExactRational& b, const ExactLimits& l,
                                               core::ResourcePlan& p) const noexcept {
    return binary_rational(*this, b, Arithmetic::Add, l, p);
}
core::Result<ExactRational> ExactRational::subtract(const ExactRational& b, const ExactLimits& l,
                                                    core::ResourcePlan& p) const noexcept {
    return binary_rational(*this, b, Arithmetic::Subtract, l, p);
}
core::Result<ExactRational> ExactRational::multiply(const ExactRational& b, const ExactLimits& l,
                                                    core::ResourcePlan& p) const noexcept {
    return binary_rational(*this, b, Arithmetic::Multiply, l, p);
}
core::Result<ExactRational> ExactRational::divide(const ExactRational& b, const ExactLimits& l,
                                                  core::ResourcePlan& p) const noexcept {
    return binary_rational(*this, b, Arithmetic::Divide, l, p);
}
core::Result<ExactRational> ExactRational::clone(const ExactLimits& limits,
                                                 core::ResourcePlan& plan) const noexcept {
    return entry(limits, [&]() -> core::Result<ExactRational> {
        const auto& a = ExactFactory::value(*this);
        check_value(a, limits);
        const auto info = admit(limits, value_bits(a), 0, retained_bytes(a), 2, Scratch::Copy);
        return ExactFactory::rational(Rational{a}, limits, plan, info);
    });
}
core::Result<std::strong_ordering> ExactRational::compare(const ExactRational& other,
                                                          const ExactLimits& limits,
                                                          core::ResourcePlan&) const noexcept {
    return entry(limits, [&]() -> core::Result<std::strong_ordering> {
        const auto& a = ExactFactory::value(*this);
        const auto& b = ExactFactory::value(other);
        check_value(a, limits);
        check_value(b, limits);
        (void)admit(limits,
                    std::max(plus(bits(a.backend().num()), bits(b.backend().denom())),
                             plus(bits(b.backend().num()), bits(a.backend().denom()))),
                    0, plus(retained_bytes(a), retained_bytes(b)), 2, Scratch::Compare);
        const auto result = a.backend().compare(b.backend());
        return result < 0   ? std::strong_ordering::less
               : result > 0 ? std::strong_ordering::greater
                            : std::strong_ordering::equal;
    });
}

namespace detail {
core::Result<PolynomialValue> make_polynomial(std::span<const ExactInputTerm> input,
                                              const ExactLimits& limits,
                                              core::ResourcePlan& plan) noexcept {
    return entry(limits, [&]() -> core::Result<PolynomialValue> {
        std::uint64_t maximum = 1, borrowed = 0;
        for (const auto& term : input) {
            check_powers(term.powers, limits);
            const auto& v = ExactFactory::value(term.coefficient);
            check_value(v, limits);
            maximum = std::max(maximum, value_bits(v));
            borrowed = plus(borrowed, retained_bytes(v));
        }
        const auto info =
            poly_plan(limits, merge_bound(input.size(), maximum), input.size(), borrowed);
        auto data = raw(input.size(), plan);
        for (std::size_t i = 0; i < input.size(); ++i) {
            data.values()[i].powers = input[i].powers;
            data.values()[i].value = ExactFactory::value(input[i].coefficient);
        }
        return finish(std::move(data), limits, plan, info);
    });
}
core::Result<PolynomialValue> clone_polynomial(const PolynomialValue& value,
                                               const ExactLimits& limits,
                                               core::ResourcePlan& plan) noexcept {
    return entry(limits, [&]() -> core::Result<PolynomialValue> {
        const auto terms = value.terms();
        const auto maximum = maximum_bits(terms, limits);
        const auto info = poly_plan(limits, merge_bound(terms.size(), maximum), terms.size(),
                                    retained_bytes(terms));
        auto data = raw(terms.size(), plan);
        for (std::size_t i = 0; i < terms.size(); ++i) {
            data.values()[i].powers = terms[i].powers;
            data.values()[i].value = ExactFactory::value(terms[i].coefficient);
        }
        return finish(std::move(data), limits, plan, info);
    });
}
core::Result<PolynomialValue> add_polynomial(const PolynomialValue& left,
                                             const PolynomialValue& right, bool subtract,
                                             const ExactLimits& limits,
                                             core::ResourcePlan& plan) noexcept {
    return entry(limits, [&]() -> core::Result<PolynomialValue> {
        const auto a = left.terms(), b = right.terms();
        const auto maximum = std::max(maximum_bits(a, limits), maximum_bits(b, limits));
        const auto count = plus(a.size(), b.size());
        const auto info = poly_plan(limits, merge_bound(count, maximum), count,
                                    plus(retained_bytes(a), retained_bytes(b)));
        auto data = raw(count, plan);
        std::size_t i = 0;
        for (const auto& t : a) {
            data.values()[i].powers = t.powers;
            data.values()[i++].value = ExactFactory::value(t.coefficient);
        }
        for (const auto& t : b) {
            data.values()[i].powers = t.powers;
            data.values()[i].value = ExactFactory::value(t.coefficient);
            if (subtract)
                data.values()[i].value.backend().negate();
            ++i;
        }
        return finish(std::move(data), limits, plan, info);
    });
}
core::Result<PolynomialValue> multiply_polynomial(const PolynomialValue& left,
                                                  const PolynomialValue& right,
                                                  const ExactLimits& limits,
                                                  core::ResourcePlan& plan) noexcept {
    return entry(limits, [&]() -> core::Result<PolynomialValue> {
        const auto a = left.terms(), b = right.terms();
        const auto maximum = plus(maximum_bits(a, limits), maximum_bits(b, limits));
        const auto count = times(a.size(), b.size());
        for (const auto& u : a)
            for (const auto& v : b)
                check_powers({index(plus(static_cast<std::uint64_t>(u.powers.first),
                                         static_cast<std::uint64_t>(v.powers.first))),
                              index(plus(static_cast<std::uint64_t>(u.powers.second),
                                         static_cast<std::uint64_t>(v.powers.second)))},
                             limits);
        const auto info = poly_plan(limits, merge_bound(count, maximum), count,
                                    plus(retained_bytes(a), retained_bytes(b)));
        auto data = raw(count, plan);
        std::size_t i = 0;
        for (const auto& u : a)
            for (const auto& v : b) {
                data.values()[i].powers = {u.powers.first + v.powers.first,
                                           u.powers.second + v.powers.second};
                data.values()[i++].value =
                    calculate(ExactFactory::value(u.coefficient),
                              ExactFactory::value(v.coefficient), Arithmetic::Multiply);
            }
        return finish(std::move(data), limits, plan, info);
    });
}
core::Result<PolynomialValue> scale_polynomial(const PolynomialValue& value,
                                               const ExactRational& scalar,
                                               const ExactLimits& limits,
                                               core::ResourcePlan& plan) noexcept {
    return entry(limits, [&]() -> core::Result<PolynomialValue> {
        const auto terms = value.terms();
        const auto& factor = ExactFactory::value(scalar);
        check_value(factor, limits);
        const auto maximum = plus(maximum_bits(terms, limits), value_bits(factor));
        const auto info = poly_plan(limits, merge_bound(terms.size(), maximum), terms.size(),
                                    plus(retained_bytes(terms), retained_bytes(factor)));
        auto data = raw(terms.size(), plan);
        for (std::size_t i = 0; i < terms.size(); ++i) {
            data.values()[i].powers = terms[i].powers;
            data.values()[i].value =
                calculate(ExactFactory::value(terms[i].coefficient), factor, Arithmetic::Multiply);
        }
        return finish(std::move(data), limits, plan, info);
    });
}
core::Result<PolynomialValue> differentiate_polynomial(const PolynomialValue& value, ExactAxis axis,
                                                       core::Index order, const ExactLimits& limits,
                                                       core::ResourcePlan& plan) noexcept {
    return entry(limits, [&]() -> core::Result<PolynomialValue> {
        check_axis(axis);
        if (order < 0)
            stop(core::FailureCode::InvalidInput);
        const auto terms = value.terms();
        const auto initial = maximum_bits(terms, limits);
        std::uint64_t degree = 0, count = 0;
        for (const auto& term : terms) {
            const auto d = axis == ExactAxis::First ? term.powers.first : term.powers.second;
            degree = std::max(degree, static_cast<std::uint64_t>(d));
            if (d >= order)
                ++count;
        }
        const auto growth =
            count == 0 ? 0 : times(static_cast<std::uint64_t>(order), std::bit_width(degree));
        const auto info = poly_plan(limits, merge_bound(count, plus(initial, growth)), count,
                                    retained_bytes(terms));
        auto data = raw(count, plan);
        std::size_t i = 0;
        for (const auto& term : terms) {
            const auto d = axis == ExactAxis::First ? term.powers.first : term.powers.second;
            if (d < order)
                continue;
            auto powers = term.powers;
            if (axis == ExactAxis::First)
                powers.first -= order;
            else
                powers.second -= order;
            Rational coefficient = ExactFactory::value(term.coefficient);
            for (core::Index k = 0; k < order; ++k)
                coefficient = calculate(coefficient, Rational{d - k}, Arithmetic::Multiply);
            data.values()[i].powers = powers;
            data.values()[i++].value = std::move(coefficient);
        }
        return finish(std::move(data), limits, plan, info);
    });
}
core::Result<PolynomialValue> integrate_polynomial(const PolynomialValue& value, ExactAxis axis,
                                                   const ExactRational& lower,
                                                   const ExactLimits& limits,
                                                   core::ResourcePlan& plan) noexcept {
    return entry(limits, [&]() -> core::Result<PolynomialValue> {
        check_axis(axis);
        const auto terms = value.terms();
        const auto& point = ExactFactory::value(lower);
        check_value(point, limits);
        const auto initial = maximum_bits(terms, limits);
        std::uint64_t degree = 0;
        for (const auto& term : terms) {
            auto powers = term.powers;
            const auto d = plus(
                static_cast<std::uint64_t>(axis == ExactAxis::First ? powers.first : powers.second),
                1);
            if (axis == ExactAxis::First)
                powers.first = index(d);
            else
                powers.second = index(d);
            check_powers(powers, limits);
            degree = std::max(degree, d);
        }
        const auto count = times(2, terms.size());
        const auto coefficient =
            plus(plus(initial, std::bit_width(degree)), times(degree, value_bits(point)));
        const auto info = poly_plan(limits, merge_bound(count, coefficient), count,
                                    plus(retained_bytes(terms), retained_bytes(point)));
        auto data = raw(count, plan);
        std::size_t i = 0;
        for (const auto& term : terms) {
            auto powers = term.powers;
            auto constant = powers;
            const auto d = axis == ExactAxis::First ? ++powers.first : ++powers.second;
            Rational coefficient_value =
                calculate(ExactFactory::value(term.coefficient), Rational{d}, Arithmetic::Divide);
            Rational constant_value =
                calculate(coefficient_value, rational_power(point, static_cast<std::uint64_t>(d)),
                          Arithmetic::Multiply);
            constant_value.backend().negate();
            if (axis == ExactAxis::First)
                constant.first = 0;
            else
                constant.second = 0;
            data.values()[i].powers = powers;
            data.values()[i++].value = std::move(coefficient_value);
            data.values()[i].powers = constant;
            data.values()[i++].value = std::move(constant_value);
        }
        return finish(std::move(data), limits, plan, info);
    });
}
core::Result<ExactRational> evaluate_polynomial(const PolynomialValue& value,
                                                const ExactRational& first,
                                                const ExactRational& second,
                                                const ExactLimits& limits,
                                                core::ResourcePlan& plan) noexcept {
    return entry(limits, [&]() -> core::Result<ExactRational> {
        const auto terms = value.terms();
        const auto& x = ExactFactory::value(first);
        const auto& y = ExactFactory::value(second);
        check_value(x, limits);
        check_value(y, limits);
        const auto initial = maximum_bits(terms, limits);
        std::uint64_t coefficient = initial;
        for (const auto& term : terms)
            coefficient = std::max(
                coefficient,
                plus(plus(initial,
                          times(static_cast<std::uint64_t>(term.powers.first), value_bits(x))),
                     times(static_cast<std::uint64_t>(term.powers.second), value_bits(y))));
        const auto info =
            poly_plan(limits, merge_bound(terms.size(), coefficient), terms.size(),
                      plus(retained_bytes(terms), plus(retained_bytes(x), retained_bytes(y))));
        Rational result = 0;
        for (const auto& term : terms) {
            const auto sx = rational_power(x, static_cast<std::uint64_t>(term.powers.first));
            const auto sy = rational_power(y, static_cast<std::uint64_t>(term.powers.second));
            const auto xy = calculate(sx, sy, Arithmetic::Multiply);
            const auto contribution =
                calculate(ExactFactory::value(term.coefficient), xy, Arithmetic::Multiply);
            result = calculate(result, contribution, Arithmetic::Add);
        }
        return ExactFactory::rational(std::move(result), limits, plan, info);
    });
}
namespace {
std::uint64_t count_power(std::uint64_t base, std::uint64_t exponent) {
    std::uint64_t result = 1;
    while (exponent != 0) {
        if ((exponent & 1) != 0)
            result = times(result, base);
        exponent >>= 1;
        if (exponent != 0)
            base = times(base, base);
    }
    return result;
}
Raw multiply_raw(const Raw& left, std::span<const ExactTerm> right, const ExactLimits& limits,
                 core::ResourcePlan& plan) {
    auto result = raw(times(left.size(), right.size()), plan);
    std::size_t i = 0;
    for (const auto& a : left.values())
        for (const auto& b : right) {
            const auto powers =
                ExactPowers{index(plus(static_cast<std::uint64_t>(a.powers.first),
                                       static_cast<std::uint64_t>(b.powers.first))),
                            index(plus(static_cast<std::uint64_t>(a.powers.second),
                                       static_cast<std::uint64_t>(b.powers.second)))};
            check_powers(powers, limits);
            result.values()[i].powers = powers;
            result.values()[i++].value =
                calculate(a.value, ExactFactory::value(b.coefficient), Arithmetic::Multiply);
        }
    return result;
}
} // namespace
core::Result<PolynomialValue> substitute_polynomial(const PolynomialValue& value,
                                                    const PolynomialValue& first,
                                                    const PolynomialValue& second,
                                                    const ExactLimits& limits,
                                                    core::ResourcePlan& plan) noexcept {
    return entry(limits, [&]() -> core::Result<PolynomialValue> {
        const auto terms = value.terms(), x = first.terms(), y = second.terms();
        const auto initial = maximum_bits(terms, limits), xb = maximum_bits(x, limits),
                   yb = maximum_bits(y, limits);
        const auto degrees = [](std::span<const ExactTerm> values) {
            std::array<std::uint64_t, 3> result{};
            for (const auto& term : values) {
                const auto a = static_cast<std::uint64_t>(term.powers.first);
                const auto b = static_cast<std::uint64_t>(term.powers.second);
                result[0] = std::max(result[0], a);
                result[1] = std::max(result[1], b);
                result[2] = std::max(result[2], plus(a, b));
            }
            return result;
        };
        const auto dx = degrees(x), dy = degrees(y);
        std::uint64_t count = 0, coefficient = initial;
        for (const auto& term : terms) {
            const auto a = static_cast<std::uint64_t>(term.powers.first),
                       b = static_cast<std::uint64_t>(term.powers.second);
            // Zero replacements annihilate this whole monomial. Validate every
            // input first, but do not build irrelevant powers or scratch plans.
            if ((x.empty() && a != 0) || (y.empty() && b != 0))
                continue;
            const auto first_degree = plus(times(a, dx[0]), times(b, dy[0]));
            const auto second_degree = plus(times(a, dx[1]), times(b, dy[1]));
            const auto total_degree = plus(times(a, dx[2]), times(b, dy[2]));
            if (first_degree > limits.maximum_axis_degree ||
                second_degree > limits.maximum_axis_degree ||
                total_degree > limits.maximum_total_degree)
                stop(core::FailureCode::ResourceExceeded);
            const auto local = times(count_power(x.size(), a), count_power(y.size(), b));
            count = plus(count, local);
            coefficient = std::max(coefficient, plus(plus(initial, times(a, xb)), times(b, yb)));
        }
        const auto info =
            poly_plan(limits, merge_bound(count, coefficient), count,
                      plus(retained_bytes(terms), plus(retained_bytes(x), retained_bytes(y))));
        auto output = raw(count, plan);
        std::size_t next = 0;
        for (const auto& term : terms) {
            if ((x.empty() && term.powers.first != 0) || (y.empty() && term.powers.second != 0))
                continue;
            auto current = raw(1, plan);
            current.values()[0].value = ExactFactory::value(term.coefficient);
            for (core::Index i = 0; i < term.powers.first; ++i)
                current = multiply_raw(current, x, limits, plan);
            for (core::Index i = 0; i < term.powers.second; ++i)
                current = multiply_raw(current, y, limits, plan);
            for (auto& result : current.values())
                output.values()[next++] = std::move(result);
        }
        return finish(std::move(output), limits, plan, info);
    });
}
} // namespace detail
} // namespace uvlas::structure
