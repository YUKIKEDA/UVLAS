#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <string_view>
#include <uvlas/core/failure.hpp>

namespace uvlas::io::detail {

// Represents the exact significand and exponent without a floating-point conversion.
class NastranRealField {
public:
    [[nodiscard]] std::string_view decimal() const& noexcept {
        return {bytes_.data(), size_};
    }
    [[nodiscard]] std::string_view decimal() const&& = delete;

private:
    friend core::Result<NastranRealField> read_nastran_real_field(std::string_view) noexcept;
    NastranRealField() = default;

    // One physical line plus a possible leading/trailing zero and implicit exponent marker.
    std::array<char, 82> bytes_{};
    std::size_t size_{};
};

// A classified Bulk Data field, bounded by the 80-byte physical line. Outer spaces
// are padding. A Real has a decimal point and may have E, D, or a signed implicit
// exponent. The result spells the exact value for the existing SI decimal handling.
[[nodiscard]] core::Result<NastranRealField>
read_nastran_real_field(std::string_view field) noexcept;

[[nodiscard]] core::Result<std::int64_t>
read_nastran_integer_field(std::string_view field) noexcept;

} // namespace uvlas::io::detail
