#pragma once

#include <concepts>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string_view>
#include <type_traits>
#include <uvlas/core/checked_dimensions.hpp>
#include <uvlas/core/failure.hpp>
#include <uvlas/core/identifier.hpp>
#include <variant>

namespace uvlas::io {

/// Why a decoded identifier or native integer cannot be converted.
enum class ScalarConversionReason : std::uint8_t {
    EmptyId,            ///< An identifier must contain at least one ASCII digit.
    InvalidIdCharacter, ///< A decoded byte is not an ASCII digit.
    IdOverflow,         ///< The next decimal digit would exceed UINT64_MAX.
    NegativeIndex,      ///< A signed native integer is below zero.
    IndexOverflow,      ///< An unsigned native integer exceeds INT64_MAX.
    UnexpectedException ///< A standard-library exception was translated at the entry boundary.
};

/// Owned, allocation-free conversion detail, independent of the source buffer.
/// Positions refer to decoded identifier bytes, never original JSON/file offsets.
struct ScalarConversionDiagnostic {
    /// Classification supplementing the returned InvalidInput Failure.
    ScalarConversionReason reason;
    /// First failing decoded byte; absent for an empty identifier or invalid index.
    std::optional<std::size_t> decoded_byte_position;
    /// Preserve the supplied reason and optional decoded position without allocation.
    ScalarConversionDiagnostic(ScalarConversionReason cause,
                               std::optional<std::size_t> position) noexcept
        : reason(cause), decoded_byte_position(position) {}
    /// Compare the owned reason and whether a decoded position is available.
    bool operator==(const ScalarConversionDiagnostic&) const = default;
};

/// An explicitly supplied native signed/unsigned 64-bit integer, not a JSON node.
/// Boolean, floating-point, string and narrower integer inputs cannot construct it.
/// A document adapter must check the node kind before constructing this value.
class NativeInteger {
public:
    /// Preserve exactly an int64_t or uint64_t, including its signedness.
    /// Construction is explicit; callers must not round a decimal node into this type.
    template <class T>
        requires(std::same_as<std::remove_cvref_t<T>, std::int64_t> ||
                 std::same_as<std::remove_cvref_t<T>, std::uint64_t>)
    explicit NativeInteger(T integer) noexcept : value_(integer) {}
    /// Return the selected native integer by value without borrowing storage.
    [[nodiscard]] std::variant<std::int64_t, std::uint64_t> selection() const noexcept {
        return value_;
    }

private:
    std::variant<std::int64_t, std::uint64_t> value_;
};

/// Existing physical identity domains accepted by decoded decimal identifier input.
/// Generation counters, accepted/trial states and array indices are excluded.
template <class T>
concept PhysicalInputId =
    std::same_as<T, core::DofId> || std::same_as<T, core::PanelId> ||
    std::same_as<T, core::WakeVertexId> || std::same_as<T, core::WakeFaceId> ||
    std::same_as<T, core::LayerId> || std::same_as<T, core::ObservationId> ||
    std::same_as<T, core::ScopeId> || std::same_as<T, core::ComponentId>;

namespace detail {
[[nodiscard]] core::Result<std::uint64_t>
read_id_digits(std::string_view text, std::optional<ScalarConversionDiagnostic>& diagnostic,
               core::FailureDetails details) noexcept;
}

/// Convert nonempty decoded ASCII digits to the requested physical identity domain.
/// Leading zeros and every uint64 value, including zero, are accepted exactly.
/// Signs, whitespace, non-ASCII bytes, embedded NULs and overflow return InvalidInput
/// with Operation::PrepareModel and the supplied context. The first failing byte wins.
/// The diagnostic is reset on entry, empty on success, and owns detail on failure.
/// It must not alias text or be shared by concurrent/reentrant calls. The input is
/// borrowed only during this call; no source/model membership is certified.
/// No allocation or exception occurs, and the input/context are not changed.
template <PhysicalInputId Id>
[[nodiscard]] core::Result<Id>
read_physical_id(std::string_view text, std::optional<ScalarConversionDiagnostic>& diagnostic,
                 core::FailureDetails details = {}) noexcept {
    const auto parsed = detail::read_id_digits(text, diagnostic, details);
    if (!parsed.has_value())
        return std::unexpected(parsed.error());
    return Id{*parsed};
}

/// Convert a native integer to an index in [0, INT64_MAX], preserving its exact value.
/// Negative values and unsigned overflow return InvalidInput/PrepareModel and an
/// owned diagnostic without a decoded byte position. Supplied FailureDetails survive.
/// The diagnostic is reset on entry and empty on success; use a distinct, nonaliased
/// slot per concurrent/reentrant call. No allocation is performed and no exception escapes.
/// This does not check an array length, validate a JSON node or retain its source.
[[nodiscard]] core::Result<core::Index>
read_index(NativeInteger input, std::optional<ScalarConversionDiagnostic>& diagnostic,
           core::FailureDetails details = {}) noexcept;

} // namespace uvlas::io
