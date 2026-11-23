#pragma once

#include "si_exact_words.hpp"

#include <uvlas/core/error_certificate.hpp>

namespace uvlas::io::detail {

/// Exact SI scale and nonnegative dimensionless coefficients for E <= a*s + r*abs(value).
struct SiRationalProjectionBudget {
    SiRationalWordsView scale, absolute, relative;
};

struct SiRationalProjectionLimits {
    std::uint64_t input_integer_bits = 0; ///< Maximum numerator or denominator extent.
    std::uint64_t steps = 0;              ///< Cap for all word inspections and arithmetic visits.
};

struct SiRationalProjection {
    SiEnclosure value;
    std::uint64_t used_steps = 0;
};

/// Project the canonical SI value itself. Recheck the retained words and outward radius.
[[nodiscard]] core::Result<SiRationalProjection>
project_rational_words(SiRationalWordsView, const SiRationalProjectionBudget&, core::Precision,
                       SiRationalProjectionLimits, const std::stop_token&,
                       core::FailureDetails) noexcept;

} // namespace uvlas::io::detail
