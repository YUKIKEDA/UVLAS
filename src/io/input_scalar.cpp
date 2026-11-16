#include <limits>
#include <uvlas/io/input_scalar.hpp>

namespace uvlas::io {
namespace {
[[nodiscard]] std::unexpected<core::Failure>
reject(ScalarConversionReason reason, std::optional<std::size_t> position,
       std::optional<ScalarConversionDiagnostic>& diagnostic,
       core::FailureDetails details) noexcept {
    diagnostic.emplace(reason, position);
    return std::unexpected(
        core::Failure{core::FailureCode::InvalidInput, core::Operation::PrepareModel, details});
}
} // namespace

core::Result<std::uint64_t>
detail::read_id_digits(std::string_view text, std::optional<ScalarConversionDiagnostic>& diagnostic,
                       core::FailureDetails details) noexcept {
    diagnostic.reset();
    if (text.empty())
        return reject(ScalarConversionReason::EmptyId, std::nullopt, diagnostic, details);
    std::uint64_t value = 0;
    constexpr auto maximum = std::numeric_limits<std::uint64_t>::max();
    for (std::size_t i = 0; i < text.size(); ++i) {
        const auto character = text[i];
        // Check bytes directly: locale digit predicates can accept another grammar.
        if (character < '0' || character > '9')
            return reject(ScalarConversionReason::InvalidIdCharacter, i, diagnostic, details);
        const auto digit = static_cast<std::uint64_t>(character - '0');
        // Subtract the known digit before multiplying: even MAX+1 never wraps.
        if (value > (maximum - digit) / 10)
            return reject(ScalarConversionReason::IdOverflow, i, diagnostic, details);
        value = value * 10 + digit;
    }
    return value;
}

core::Result<core::Index> read_index(NativeInteger input,
                                     std::optional<ScalarConversionDiagnostic>& diagnostic,
                                     core::FailureDetails details) noexcept {
    diagnostic.reset();
    const auto selected = input.selection();
    const auto result =
        core::guard_entry(core::Operation::PrepareModel, core::FailureCode::InvalidInput, details,
                          [&]() -> core::Result<core::Index> {
                              // Primitive alternatives cannot become valueless. Still contain the
                              // standard visitor exception path at the library boundary.
                              return std::visit(
                                  [&](auto value) -> core::Result<core::Index> {
                                      using Integer = decltype(value);
                                      if constexpr (std::same_as<Integer, std::int64_t>) {
                                          // Reject before any unsigned conversion or negation,
                                          // including INT64_MIN.
                                          if (value < 0)
                                              return reject(ScalarConversionReason::NegativeIndex,
                                                            std::nullopt, diagnostic, details);
                                          return value;
                                      } else {
                                          constexpr auto maximum = static_cast<std::uint64_t>(
                                              std::numeric_limits<core::Index>::max());
                                          if (value > maximum)
                                              return reject(ScalarConversionReason::IndexOverflow,
                                                            std::nullopt, diagnostic, details);
                                          // The range check makes this narrowing conversion exact.
                                          return static_cast<core::Index>(value);
                                      }
                                  },
                                  selected);
                          });
    if (!result.has_value() && !diagnostic.has_value())
        diagnostic.emplace(ScalarConversionReason::UnexpectedException, std::nullopt);
    return result;
}
} // namespace uvlas::io
