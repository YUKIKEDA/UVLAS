#pragma once

#include <compare>
#include <cstdint>

namespace uvlas::core {

// All uint64 values are legal. Absence is an optional, never a reserved ID.
template <class Tag>
class Identifier {
public:
    explicit constexpr Identifier(std::uint64_t value) noexcept : value_(value) {}
    [[nodiscard]] constexpr std::uint64_t value() const noexcept {
        return value_;
    }
    auto operator<=>(const Identifier&) const = default;

private:
    std::uint64_t value_;
};

namespace detail {
struct DofTag;
struct PanelTag;
struct WakeVertexTag;
struct WakeFaceTag;
struct LayerTag;
struct TrialTag;
struct ModelGenerationTag;
struct NumericalGenerationTag;
struct AcceptedStateTag;
struct ObservationTag;
struct ScopeTag;
struct ComponentTag;
} // namespace detail

using DofId = Identifier<detail::DofTag>;
using PanelId = Identifier<detail::PanelTag>;
using WakeVertexId = Identifier<detail::WakeVertexTag>;
using WakeFaceId = Identifier<detail::WakeFaceTag>;
using LayerId = Identifier<detail::LayerTag>;
using TrialId = Identifier<detail::TrialTag>;
using ModelGeneration = Identifier<detail::ModelGenerationTag>;
using NumericalGeneration = Identifier<detail::NumericalGenerationTag>;
using AcceptedStateId = Identifier<detail::AcceptedStateTag>;
using ObservationId = Identifier<detail::ObservationTag>;
using ScopeId = Identifier<detail::ScopeTag>;
using ComponentId = Identifier<detail::ComponentTag>;

} // namespace uvlas::core
