#pragma once

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/rational_adaptor.hpp>

#include <cstddef>
#include <uvlas/structure/exact_polynomial.hpp>

namespace uvlas::structure::detail {
using IntegerBackend = boost::multiprecision::cpp_int_backend<>;
using Integer = boost::multiprecision::number<IntegerBackend, boost::multiprecision::et_off>;
using Rational =
    boost::multiprecision::number<boost::multiprecision::rational_adaptor<IntegerBackend>,
                                  boost::multiprecision::et_off>;
/// Payload layout used by callers that preflight cumulative exact-algebra arrays.
[[nodiscard]] std::size_t exact_polynomial_raw_term_bytes() noexcept;
/// Check the existing rational-operation storage envelope without materializing any operands.
[[nodiscard]] core::Result<ExactPreparationInfo>
exact_rational_admission(const ExactLimits&, std::uint64_t work_bits,
                         std::uint64_t borrowed_bytes) noexcept;
struct RationalStorage {
    Rational value;
};
struct PolynomialStorage {
    core::PlannedBuffer<ExactTerm> terms;
};
struct ExactFactory {
    static const Rational& value(const ExactRational&);
    static core::Result<ExactRational> rational(Rational, const ExactLimits&, core::ResourcePlan&,
                                                ExactPreparationInfo);
    static core::Result<PolynomialValue> polynomial(core::PlannedBuffer<ExactTerm>,
                                                    core::ResourcePlan&, ExactPreparationInfo);
    static void timing(ExactRational&, std::optional<core::StageTiming>) noexcept;
    static void timing(PolynomialValue&, std::optional<core::StageTiming>) noexcept;
};
} // namespace uvlas::structure::detail
