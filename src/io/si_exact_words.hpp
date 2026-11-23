#pragma once

#include <span>
#include <uvlas/io/si_case.hpp>

namespace uvlas::io::detail {

struct SiRationalWordsView {
    int sign;
    std::span<const std::uint32_t> numerator, denominator;
};

/// Canonical little-endian 32-bit magnitudes for an exact SI recipe.
class SiRationalWords {
public:
    SiRationalWords() noexcept = default;
    SiRationalWords(SiRationalWords&&) noexcept = default;
    SiRationalWords& operator=(SiRationalWords&&) noexcept = default;
    SiRationalWords(const SiRationalWords&) = delete;
    SiRationalWords& operator=(const SiRationalWords&) = delete;

    [[nodiscard]] SiRationalWordsView view() const& noexcept;
    [[nodiscard]] SiRationalWordsView view() const&& = delete;

private:
    friend core::Result<SiRationalWords>
    materialize_rational_words(const SiRecipeView&, const structure::ExactLimits&,
                               core::ResourcePlan&, const std::stop_token&) noexcept;
    SiRationalWords(core::PlannedBuffer<std::uint32_t> words, std::size_t numerator,
                    int sign) noexcept;
    core::PlannedBuffer<std::uint32_t> words_;
    std::size_t numerator_ = 0;
    int sign_ = 0;
};

[[nodiscard]] core::Result<SiRationalWords>
materialize_rational_words(const SiRecipeView&, const structure::ExactLimits&, core::ResourcePlan&,
                           const std::stop_token&) noexcept;

} // namespace uvlas::io::detail
