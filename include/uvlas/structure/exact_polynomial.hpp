#pragma once

#include <compare>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string_view>
#include <utility>
#include <uvlas/core/instrumentation.hpp>
#include <uvlas/core/resource_plan.hpp>
#include <variant>

namespace uvlas::structure {
/// Admit a source-derived estimate of simultaneous external integer-limb payload.
/// The cap is not an allocator limit or a process-resident-memory guarantee.
struct ExactEstimateAdmission {
    std::uint64_t bytes = 0; ///< Explicit cap in bytes; zero is not unlimited;
};
/// Require a proven external allocation bound; currently returns Unsupported.
struct ExactRequireProvenBound {
    std::uint64_t bytes = 0; ///< Explicit cap in bytes; zero is not unlimited;
};
/// Select the external-memory admission policy explicitly for each preparation.
using ExactMemoryLimit = std::variant<ExactEstimateAdmission, ExactRequireProvenBound>;
/// Enable product-owned preparation timing without changing exact numerical values.
enum class ExactTimingMode : std::uint8_t {
    Disabled, ///< Read no clock and omit timing diagnostics.
    Enabled   ///< Measure the producing preparation call.
};
/// Explicit limits; zero never means unlimited. Input/output values and temporary
/// integers have distinct bit caps. Polynomial raw terms are checked before merging.
struct ExactLimits {
    /// Configure every cap explicitly; values are validated at each operation.
    ExactLimits(std::uint64_t input_characters, std::uint64_t terms, std::uint64_t raw_terms,
                std::uint64_t axis_degree, std::uint64_t total_degree, std::uint64_t value_bits,
                std::uint64_t work_bits, ExactMemoryLimit memory,
                ExactTimingMode measurement) noexcept
        : maximum_input_characters(input_characters),
          maximum_terms(terms),
          maximum_raw_terms(raw_terms),
          maximum_axis_degree(axis_degree),
          maximum_total_degree(total_degree),
          maximum_value_bits(value_bits),
          maximum_work_bits(work_bits),
          external_memory(memory),
          timing(measurement) {}
    std::uint64_t maximum_input_characters; ///< Maximum bytes in one numeric input string.
    std::uint64_t maximum_terms;            ///< Maximum nonzero canonical polynomial terms.
    std::uint64_t maximum_raw_terms;        ///< Maximum raw terms before exact merging.
    std::uint64_t maximum_axis_degree;      ///< Maximum degree in either coordinate.
    std::uint64_t maximum_total_degree;     ///< Maximum sum of the two powers, within Index.
    std::uint64_t maximum_value_bits;       ///< Maximum numerator/denominator bits of owned values.
    std::uint64_t maximum_work_bits;  ///< Guarded, word-rounded temporary-integer storage bit cap.
    ExactMemoryLimit external_memory; ///< Limb-payload admission, distinct from ResourcePlan.
    ExactTimingMode timing;           ///< Optional product preparation timing.
};
/// Diagnostics for the producing call, independent of its ResourcePlan lifetime.
/// External bytes and work bits are forecasts, not observed allocations or RSS.
struct ExactPreparationInfo {
    std::uint64_t estimated_external_payload_bytes = 0;
    std::uint64_t algebraic_work_bits = 0; ///< Bound before backend word rounding and guards.
    std::uint64_t planned_work_bits = 0;   ///< Guarded storage envelope; may differ by platform.
    std::uint64_t planned_raw_terms = 0;
    std::optional<core::StageTiming> timing;
};
namespace detail {
struct RationalStorage;
struct PolynomialStorage;
struct ExactFactory;
} // namespace detail

/// Owned reduced rational p/q with positive denominator and canonical zero 0/1.
/// This is reference-preparation algebra, not a rounded physical-state buffer.
/// Moves transfer storage and leave an exact zero; borrowed views expire on moves.
class ExactRational {
public:
    /// Construct the exact zero without allocation.
    ExactRational() noexcept;
    /// Transfer ownership, leaving the source equal to zero.
    ExactRational(ExactRational&&) noexcept;
    /// Replace this value, invalidating its old borrowed views.
    ExactRational& operator=(ExactRational&&) noexcept;
    /// Implicit copying cannot bypass the explicit allocation budget.
    ExactRational(const ExactRational&) = delete;
    ExactRational& operator=(const ExactRational&) = delete;
    /// Release owned value storage; preparation payload is not refunded.
    ~ExactRational();

    /// Parse two signed base-10 integer strings. Leading zeros and signs are
    /// accepted; whitespace, exponents, empty strings and denominator zero fail.
    [[nodiscard]] static core::Result<ExactRational>
    from_integer_pair(std::string_view numerator, std::string_view denominator, const ExactLimits&,
                      core::ResourcePlan&) noexcept;
    /// Parse a finite decimal: optional sign, decimal digits with optional dot,
    /// and optional e/E exponent. No whitespace or nonfinite spelling is accepted.
    [[nodiscard]] static core::Result<ExactRational>
    from_decimal(std::string_view, const ExactLimits&, core::ResourcePlan&) noexcept;
    /// Preserve the exact IEEE binary64 word, including subnormals. Both signed
    /// zeros become 0/1; NaN and infinities return Nonfinite.
    [[nodiscard]] static core::Result<ExactRational> from_binary64(double, const ExactLimits&,
                                                                   core::ResourcePlan&) noexcept;
    /// Copy this value into new owned storage under explicit budgets.
    [[nodiscard]] core::Result<ExactRational> clone(const ExactLimits&,
                                                    core::ResourcePlan&) const noexcept;
    /// Return this plus, minus, times or divided by the other exact value.
    /// Division by zero returns InvalidInput. Neither operand is modified.
    [[nodiscard]] core::Result<ExactRational> add(const ExactRational&, const ExactLimits&,
                                                  core::ResourcePlan&) const noexcept;
    /// Subtract the other rational exactly; neither operand changes.
    [[nodiscard]] core::Result<ExactRational> subtract(const ExactRational&, const ExactLimits&,
                                                       core::ResourcePlan&) const noexcept;
    /// Multiply exact canonical values under the preparation budgets.
    [[nodiscard]] core::Result<ExactRational> multiply(const ExactRational&, const ExactLimits&,
                                                       core::ResourcePlan&) const noexcept;
    /// Divide by a nonzero rational; zero divisors return InvalidInput.
    [[nodiscard]] core::Result<ExactRational> divide(const ExactRational&, const ExactLimits&,
                                                     core::ResourcePlan&) const noexcept;
    /// Compare exact values under an external cross-product working budget.
    [[nodiscard]] core::Result<std::strong_ordering>
    compare(const ExactRational&, const ExactLimits&, core::ResourcePlan&) const noexcept;
    /// Compare canonical components without allocation or floating-point conversion.
    [[nodiscard]] bool operator==(const ExactRational&) const noexcept;
    /// Sign of the exact numerator: -1, 0 or +1.
    [[nodiscard]] int sign() const noexcept;
    /// Bit length of the absolute numerator; zero has length zero.
    [[nodiscard]] std::uint64_t numerator_bits() const noexcept;
    /// Bit length of the positive denominator, including one for canonical zero.
    [[nodiscard]] std::uint64_t denominator_bits() const noexcept;
    /// Borrow diagnostics until this value is moved, assigned or destroyed.
    [[nodiscard]] const ExactPreparationInfo& preparation() const& noexcept;
    [[nodiscard]] const ExactPreparationInfo& preparation() const&& = delete;

private:
    friend struct detail::ExactFactory;
    explicit ExactRational(core::PlannedBuffer<detail::RationalStorage>,
                           ExactPreparationInfo) noexcept;
    core::PlannedBuffer<detail::RationalStorage> storage_;
    ExactPreparationInfo info_;
};

/// Distinct coordinate frames; algebra never infers geometric scale factors.
enum class ExactCoordinates : std::uint8_t {
    ChartST, ///< Global coordinates s,t of one graph chart.
    CellUV   ///< Local coordinates u,v of one cell.
};
/// Select a polynomial's first or second coordinate explicitly.
enum class ExactAxis : std::uint8_t {
    First, ///< Differentiate or integrate the first coordinate.
    Second ///< Differentiate or integrate the second coordinate.
};
/// Nonnegative powers, validated before conversion or arithmetic on their indices.
struct ExactPowers {
    core::Index first = 0;  ///< Power of the first coordinate.
    core::Index second = 0; ///< Power of the second coordinate.
};
/// A borrowed input term. Its coefficient must outlive the synchronous factory call.
struct ExactInputTerm {
    ExactPowers powers;
    const ExactRational& coefficient;
};
/// A stored term owns its coefficient; polynomial access exposes read-only terms.
struct ExactTerm {
    ExactPowers powers{0, 0};
    ExactRational coefficient;
};
namespace detail {
class PolynomialValue {
public:
    PolynomialValue() noexcept;
    PolynomialValue(PolynomialValue&&) noexcept;
    PolynomialValue& operator=(PolynomialValue&&) noexcept;
    PolynomialValue(const PolynomialValue&) = delete;
    PolynomialValue& operator=(const PolynomialValue&) = delete;
    ~PolynomialValue();
    [[nodiscard]] std::span<const ExactTerm> terms() const& noexcept;
    [[nodiscard]] std::span<const ExactTerm> terms() const&& = delete;
    [[nodiscard]] const ExactPreparationInfo& preparation() const& noexcept;
    [[nodiscard]] const ExactPreparationInfo& preparation() const&& = delete;

private:
    friend struct ExactFactory;
    explicit PolynomialValue(core::PlannedBuffer<PolynomialStorage>, ExactPreparationInfo) noexcept;
    core::PlannedBuffer<PolynomialStorage> storage_;
    ExactPreparationInfo info_;
};
[[nodiscard]] core::Result<PolynomialValue>
make_polynomial(std::span<const ExactInputTerm>, const ExactLimits&, core::ResourcePlan&) noexcept;
[[nodiscard]] core::Result<PolynomialValue>
clone_polynomial(const PolynomialValue&, const ExactLimits&, core::ResourcePlan&) noexcept;
[[nodiscard]] core::Result<PolynomialValue> add_polynomial(const PolynomialValue&,
                                                           const PolynomialValue&, bool subtract,
                                                           const ExactLimits&,
                                                           core::ResourcePlan&) noexcept;
[[nodiscard]] core::Result<PolynomialValue> multiply_polynomial(const PolynomialValue&,
                                                                const PolynomialValue&,
                                                                const ExactLimits&,
                                                                core::ResourcePlan&) noexcept;
[[nodiscard]] core::Result<PolynomialValue> scale_polynomial(const PolynomialValue&,
                                                             const ExactRational&,
                                                             const ExactLimits&,
                                                             core::ResourcePlan&) noexcept;
[[nodiscard]] core::Result<PolynomialValue> differentiate_polynomial(const PolynomialValue&,
                                                                     ExactAxis, core::Index,
                                                                     const ExactLimits&,
                                                                     core::ResourcePlan&) noexcept;
[[nodiscard]] core::Result<PolynomialValue> integrate_polynomial(const PolynomialValue&, ExactAxis,
                                                                 const ExactRational&,
                                                                 const ExactLimits&,
                                                                 core::ResourcePlan&) noexcept;
[[nodiscard]] core::Result<ExactRational>
evaluate_polynomial(const PolynomialValue&, const ExactRational&, const ExactRational&,
                    const ExactLimits&, core::ResourcePlan&) noexcept;
[[nodiscard]] core::Result<PolynomialValue>
substitute_polynomial(const PolynomialValue&, const PolynomialValue&, const PolynomialValue&,
                      const ExactLimits&, core::ResourcePlan&) noexcept;
} // namespace detail

/// A finite exact polynomial with an explicit coordinate frame. Terms are ordered
/// by first power and then second power, duplicates are merged and zeros removed.
/// A moved-from value is the zero polynomial. No PG basis qualification is implied.
template <ExactCoordinates Coordinates>
    requires(Coordinates == ExactCoordinates::ChartST || Coordinates == ExactCoordinates::CellUV)
class ExactPolynomial2 {
public:
    /// Construct zero without allocation; its degree is absent.
    ExactPolynomial2() noexcept = default;
    /// Transfer ownership, leaving the source equal to zero.
    ExactPolynomial2(ExactPolynomial2&&) noexcept = default;
    ExactPolynomial2& operator=(ExactPolynomial2&&) noexcept = default;
    ExactPolynomial2(const ExactPolynomial2&) = delete;
    ExactPolynomial2& operator=(const ExactPolynomial2&) = delete;
    /// Normalize borrowed terms into an independent canonical owned polynomial.
    [[nodiscard]] static core::Result<ExactPolynomial2> make(std::span<const ExactInputTerm> terms,
                                                             const ExactLimits& limits,
                                                             core::ResourcePlan& plan) noexcept {
        return wrap(detail::make_polynomial(terms, limits, plan));
    }
    /// Copy this polynomial under explicit payload and external budgets.
    [[nodiscard]] core::Result<ExactPolynomial2> clone(const ExactLimits& limits,
                                                       core::ResourcePlan& plan) const noexcept {
        return wrap(detail::clone_polynomial(value_, limits, plan));
    }
    /// Add values in this same coordinate frame without modifying either operand.
    [[nodiscard]] core::Result<ExactPolynomial2> add(const ExactPolynomial2& other,
                                                     const ExactLimits& limits,
                                                     core::ResourcePlan& plan) const noexcept {
        return wrap(detail::add_polynomial(value_, other.value_, false, limits, plan));
    }
    /// Subtract values in this same coordinate frame.
    [[nodiscard]] core::Result<ExactPolynomial2> subtract(const ExactPolynomial2& other,
                                                          const ExactLimits& limits,
                                                          core::ResourcePlan& plan) const noexcept {
        return wrap(detail::add_polynomial(value_, other.value_, true, limits, plan));
    }
    /// Multiply and normalize all raw products under their separate term cap.
    [[nodiscard]] core::Result<ExactPolynomial2> multiply(const ExactPolynomial2& other,
                                                          const ExactLimits& limits,
                                                          core::ResourcePlan& plan) const noexcept {
        return wrap(detail::multiply_polynomial(value_, other.value_, limits, plan));
    }
    /// Multiply every coefficient by an exact rational scalar.
    [[nodiscard]] core::Result<ExactPolynomial2> scale(const ExactRational& scalar,
                                                       const ExactLimits& limits,
                                                       core::ResourcePlan& plan) const noexcept {
        return wrap(detail::scale_polynomial(value_, scalar, limits, plan));
    }
    /// Differentiate by the requested coordinate and nonnegative order. Order zero
    /// clones; an order greater than the degree returns exact zero.
    [[nodiscard]] core::Result<ExactPolynomial2>
    differentiate(ExactAxis axis, core::Index order, const ExactLimits& limits,
                  core::ResourcePlan& plan) const noexcept {
        return wrap(detail::differentiate_polynomial(value_, axis, order, limits, plan));
    }
    /// Integrate from an explicit lower coordinate to the corresponding variable;
    /// the other coordinate remains symbolic. No implicit integration constant.
    [[nodiscard]] core::Result<ExactPolynomial2>
    integrate(ExactAxis axis, const ExactRational& lower, const ExactLimits& limits,
              core::ResourcePlan& plan) const noexcept {
        return wrap(detail::integrate_polynomial(value_, axis, lower, limits, plan));
    }
    /// Evaluate at two exact coordinates without floating-point conversion.
    [[nodiscard]] core::Result<ExactRational> evaluate(const ExactRational& first,
                                                       const ExactRational& second,
                                                       const ExactLimits& limits,
                                                       core::ResourcePlan& plan) const noexcept {
        return detail::evaluate_polynomial(value_, first, second, limits, plan);
    }
    /// Substitute polynomials in an explicitly selected destination frame. Both
    /// replacements share that frame; cross-frame addition remains ill-formed.
    template <ExactCoordinates Destination>
    [[nodiscard]] core::Result<ExactPolynomial2<Destination>>
    substitute(const ExactPolynomial2<Destination>& first,
               const ExactPolynomial2<Destination>& second, const ExactLimits& limits,
               core::ResourcePlan& plan) const noexcept {
        return ExactPolynomial2<Destination>::wrap(
            detail::substitute_polynomial(value_, first.value_, second.value_, limits, plan));
    }
    /// Borrow immutable terms until ownership is moved, assigned or destroyed.
    [[nodiscard]] std::span<const ExactTerm> terms() const& noexcept {
        return value_.terms();
    }
    [[nodiscard]] std::span<const ExactTerm> terms() const&& = delete;
    /// Total degree, or no degree for the zero polynomial.
    [[nodiscard]] std::optional<core::Index> degree() const noexcept {
        std::optional<core::Index> result;
        for (const auto& term : value_.terms()) {
            const auto d = term.powers.first + term.powers.second;
            if (!result || d > *result)
                result = d;
        }
        return result;
    }
    /// Borrow the producing call's allocation forecast and optional timing.
    [[nodiscard]] const ExactPreparationInfo& preparation() const& noexcept {
        return value_.preparation();
    }
    [[nodiscard]] const ExactPreparationInfo& preparation() const&& = delete;

private:
    template <ExactCoordinates Other>
        requires(Other == ExactCoordinates::ChartST || Other == ExactCoordinates::CellUV)
    friend class ExactPolynomial2;
    explicit ExactPolynomial2(detail::PolynomialValue value) noexcept : value_(std::move(value)) {}
    static core::Result<ExactPolynomial2>
    wrap(core::Result<detail::PolynomialValue> value) noexcept {
        if (!value)
            return std::unexpected(value.error());
        return ExactPolynomial2{std::move(*value)};
    }
    detail::PolynomialValue value_;
};
/// Exact polynomials expressed in chart coordinates s,t.
using ExactChartPolynomial = ExactPolynomial2<ExactCoordinates::ChartST>;
/// Exact polynomials expressed in local cell coordinates u,v.
using ExactCellPolynomial = ExactPolynomial2<ExactCoordinates::CellUV>;
} // namespace uvlas::structure
