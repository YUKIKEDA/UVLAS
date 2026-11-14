#pragma once

#include <compare>
#include <cstdint>

namespace uvlas::core {

/// Strongly typed identity for a model object or state, not an array index.
/// @tparam Tag Distinguishes identity domains even when their numbers coincide.
/// All uint64 values, including zero, are legal; represent absence with optional.
/// Construction does not verify membership in a model or a generation.
template <class Tag>
class Identifier {
public:
    /// Preserve the full identity value without interpreting it as an array offset.
    explicit constexpr Identifier(std::uint64_t value) noexcept : value_(value) {}
    /// Return the stored number by value; no model-membership check is implied.
    [[nodiscard]] constexpr std::uint64_t value() const noexcept {
        return value_;
    }
    /// Compare numbers within this identity domain only.
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

/// Identity of a structural degree of freedom.
using DofId = Identifier<detail::DofTag>;
/// Identity of an aerodynamic panel.
using PanelId = Identifier<detail::PanelTag>;
/// Identity of a wake vertex.
using WakeVertexId = Identifier<detail::WakeVertexTag>;
/// Identity of a wake face.
using WakeFaceId = Identifier<detail::WakeFaceTag>;
/// Identity of a material layer.
using LayerId = Identifier<detail::LayerTag>;
/// Identity of a trial state; distinct from an accepted state.
using TrialId = Identifier<detail::TrialTag>;
/// Generation of the model definition.
using ModelGeneration = Identifier<detail::ModelGenerationTag>;
/// Generation of the numerical state used for evaluation.
using NumericalGeneration = Identifier<detail::NumericalGenerationTag>;
/// Identity of an accepted state, used for both parent and last-accepted context.
using AcceptedStateId = Identifier<detail::AcceptedStateTag>;
/// Identity of the measured quantity whose error is being assessed.
using ObservationId = Identifier<detail::ObservationTag>;
/// Identity of the domain or subset over which an error bound was evaluated.
using ScopeId = Identifier<detail::ScopeTag>;
/// Identity of the component associated with a diagnostic.
using ComponentId = Identifier<detail::ComponentTag>;

} // namespace uvlas::core
