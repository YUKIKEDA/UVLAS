#include "nastran_number_fields.hpp"

#include <algorithm>
#include <charconv>
#include <span>

namespace uvlas::io::detail {
namespace {

constexpr std::size_t physical_line_width = 80;

[[nodiscard]] auto failure(core::FailureCode code = core::FailureCode::InvalidInput) noexcept {
    return std::unexpected(core::Failure{code, core::Operation::PrepareModel});
}

[[nodiscard]] std::string_view trim(std::string_view field) noexcept {
    const auto first = field.find_first_not_of(' ');
    if (first == std::string_view::npos) {
        return {};
    }
    return field.substr(first, field.find_last_not_of(' ') - first + 1);
}

[[nodiscard]] std::string_view digits(std::string_view& field) noexcept {
    const auto end = field.find_first_not_of("0123456789");
    const auto prefix = field.substr(0, end);
    field.remove_prefix(prefix.size());
    return prefix;
}

[[nodiscard]] char sign(std::string_view& field) noexcept {
    if (field.empty() || (field.front() != '-' && field.front() != '+')) {
        return '\0';
    }
    const auto value = field.front();
    field.remove_prefix(1);
    return value;
}

} // namespace

core::Result<NastranRealField> read_nastran_real_field(std::string_view field) noexcept {
    if (field.size() > physical_line_width) {
        return failure(core::FailureCode::ResourceExceeded);
    }
    field = trim(field);
    const auto mantissa_sign = sign(field);
    const auto whole = digits(field);
    if (field.empty() || field.front() != '.') {
        return failure();
    }
    field.remove_prefix(1);
    const auto fraction = digits(field);
    if (whole.empty() && fraction.empty()) {
        return failure();
    }

    char exponent_sign = '\0';
    std::string_view exponent;
    if (!field.empty()) {
        const auto marker = field.front();
        if (marker == 'E' || marker == 'e' || marker == 'D' || marker == 'd') {
            field.remove_prefix(1);
        } else if (marker != '+' && marker != '-') {
            return failure();
        }
        exponent_sign = sign(field);
        exponent = digits(field);
        if (exponent.empty() || !field.empty()) {
            return failure();
        }
    }

    NastranRealField value;
    auto sink = std::span{value.bytes_}.begin();
    const auto put = [&](std::string_view piece) { sink = std::ranges::copy(piece, sink).out; };
    if (mantissa_sign == '-') {
        put("-");
    }
    const auto nonzero = whole.find_first_not_of('0');
    put(nonzero == std::string_view::npos ? "0" : whole.substr(nonzero));
    put(".");
    put(fraction.empty() ? "0" : fraction);
    if (!exponent.empty()) {
        put("e");
        if (exponent_sign != '\0') {
            put({&exponent_sign, 1});
        }
        put(exponent);
    }
    value.size_ = static_cast<std::size_t>(sink - std::span{value.bytes_}.begin());
    return value;
}

core::Result<std::int64_t> read_nastran_integer_field(std::string_view field) noexcept {
    if (field.size() > physical_line_width) {
        return failure(core::FailureCode::ResourceExceeded);
    }
    field = trim(field);
    auto input = field;
    if (sign(input) == '+') {
        field = input;
        if (input.starts_with('-')) {
            return failure();
        }
    }
    if (field.empty()) {
        return failure();
    }
    std::int64_t value{};
    const auto end = field.data() + field.size();
    const auto parsed = std::from_chars(field.data(), end, value);
    if (parsed.ec != std::errc{} || parsed.ptr != end) {
        return failure();
    }
    return value;
}

} // namespace uvlas::io::detail
