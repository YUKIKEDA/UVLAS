#pragma once

#include <array>
#include <optional>
#include <span>
#include <stop_token>
#include <string_view>
#include <uvlas/core/identifier.hpp>
#include <uvlas/core/resource_plan.hpp>
#include <variant>

namespace uvlas::io::detail {

// Lower and Upper refer to the selected layer's through-thickness interval.
enum class VisualizationLayerSide : std::uint8_t { Interior, Lower, Upper };
enum class VisualizationMeasure : std::uint8_t { Green, SecondPk, FirstPk, Cauchy };
enum class VisualizationBasis : std::uint8_t { World, Material };
enum class VisualizationContribution : std::uint8_t { Total, Elastic, Viscous };

struct VisualizationPartKey {
    core::ComponentId component{0};
    core::LayerId layer{0};
    VisualizationLayerSide side = VisualizationLayerSide::Interior;

    auto operator<=>(const VisualizationPartKey&) const = default;
};

struct VisualizationPoint {
    std::array<double, 3> reference_metres{};
    std::array<double, 3> current_metres{};
    std::array<double, 3> reference_normal{};
    std::array<double, 3> current_normal{};
    double reference_thickness_metres = 0;
    double current_thickness_metres = 0;
    double layer_fraction = 0;
};

using VisualizationSymmetric = std::array<double, 6>;
using VisualizationMatrix = std::array<double, 9>;
using VisualizationTriangle = std::array<std::uint64_t, 3>;

struct VisualizationTensorField {
    VisualizationMeasure measure = VisualizationMeasure::Green;
    VisualizationBasis basis = VisualizationBasis::World;
    VisualizationContribution contribution = VisualizationContribution::Total;
    // One tuple per sample point. Symmetric order is xx, yy, zz, yz, xz, xy;
    // First PK has nine row-major current-force/reference-normal components.
    std::variant<std::span<const VisualizationSymmetric>, std::span<const VisualizationMatrix>>
        values;
};

struct VisualizationPart {
    VisualizationPartKey key;
    std::span<const VisualizationPoint> points;
    std::span<const VisualizationTriangle> triangles;  // Local visualization-point indices.
    std::span<const VisualizationTensorField> tensors; // Available qualified fields only.
};

struct VisualizationFrame {
    double time_seconds = 0;
    core::AcceptedStateId accepted{0};
    std::uint64_t segment = 0;
    // Stable, ordered component/layer/side keys within the segment. Each part owns
    // a distinct logical grid even where geometrical positions coincide.
    std::span<const VisualizationPart> parts;
};

struct VisualizationLimits {
    std::size_t maximum_frames;
    std::size_t maximum_files;
    std::size_t maximum_points;
    std::size_t maximum_triangles;
    std::size_t maximum_output_bytes;
};

struct VisualizationFileView {
    std::string_view name;
    std::string_view text;
};

class VisualizationSeries {
public:
    VisualizationSeries(VisualizationSeries&&) noexcept = default;
    VisualizationSeries& operator=(VisualizationSeries&&) noexcept = default;

    [[nodiscard]] std::size_t size() const noexcept;
    [[nodiscard]] std::optional<VisualizationFileView> at(std::size_t) const& noexcept;
    [[nodiscard]] std::optional<VisualizationFileView> at(std::size_t) const&& = delete;

private:
    struct FileLocation {
        std::array<char, 64> name{};
        std::size_t name_size = 0;
        std::size_t offset = 0;
        std::size_t size = 0;
    };

    friend core::Result<VisualizationSeries>
    encode_visualization_series(std::span<const VisualizationFrame>, VisualizationLimits,
                                core::ResourcePlan&, const std::stop_token&) noexcept;

    VisualizationSeries(core::PlannedBuffer<FileLocation> files,
                        core::PlannedBuffer<char> bytes) noexcept;
    core::PlannedBuffer<FileLocation> files_;
    core::PlannedBuffer<char> bytes_;
};

// The caller supplies accepted, qualified binary64 visualization samples in world geometry
// coordinates and SI units. The encoder does not recover or project numerical fields.
// Each segment receives a PVD, followed by its independent part VTUs. Empty input
// produces no files. Failure returns no partial series.
[[nodiscard]] core::Result<VisualizationSeries>
encode_visualization_series(std::span<const VisualizationFrame>, VisualizationLimits,
                            core::ResourcePlan&, const std::stop_token& = {}) noexcept;

} // namespace uvlas::io::detail
