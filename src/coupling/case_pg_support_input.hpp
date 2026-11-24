#pragma once

#include <type_traits>
#include <uvlas/io/si_case.hpp>

namespace uvlas::coupling::detail {

/// Check the fixed reference motion spelling with at most six typed SI lookups.
[[nodiscard]] inline core::Result<bool>
stationary_case_pg_motion(const io::SiCase& source, const io::schema::Motion& motion) noexcept {
    namespace schema = io::schema;
    if (motion.events.size() != 0)
        return false;
    const auto zero_series = [&](const auto& field) -> core::Result<bool> {
        const auto series = field.record();
        const auto* constant = std::get_if<schema::TimeSeriesConstant>(&series);
        if (constant == nullptr || constant->value.size() != 3)
            return false;
        for (core::Index i = 0; i < 3; ++i) {
            const auto input = constant->value.at(i);
            if (!input)
                return std::unexpected(input.error());
            const auto quantity = source.quantity(*input);
            if (!quantity)
                return std::unexpected(quantity.error());
            const auto zero = std::visit(
                [](const auto& value) {
                    if constexpr (std::is_same_v<std::remove_cvref_t<decltype(value)>,
                                                 io::SiDecimalView>)
                        return value.sign == 0;
                    else
                        return value.numerator.sign == 0;
                },
                quantity->recipe.value);
            if (!zero)
                return false;
        }
        return true;
    };
    const auto translation = zero_series(motion.translation);
    if (!translation || !*translation)
        return translation;
    return zero_series(motion.rotation_vector);
}

} // namespace uvlas::coupling::detail
