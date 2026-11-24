#pragma once

#include "nastran_shell_cards.hpp"
#include "nastran_support_cards.hpp"

#include <uvlas/io/si_case.hpp>

namespace uvlas::io::detail {

class NastranIsotropicShellCards;

struct NastranCartesianInterval {
    structure::ExactRational lower, upper;
};

struct NastranGridEnclosure {
    NastranGridId grid{0};
    NastranCoordinateId cp{0};
    NastranBulkPosition source;
    std::optional<NastranBulkPosition> cp_field;
    std::array<std::optional<NastranBulkPosition>, 3> tuple_fields;
    std::array<NastranCartesianInterval, 3> metres;
};

enum class NastranPositionFailure : std::uint8_t {
    Source,
    Selection,
    CoordinateDefinition,
    Tuple,
    FrameOrigin,
    FrameAxis,
    FramePlane,
    Conversion,
    SpcTerm,
    DisplacementAxis,
};

struct NastranPositionDiagnostic {
    NastranPositionFailure reason = NastranPositionFailure::Source;
    std::optional<NastranGridId> grid;
    std::optional<NastranCoordinateId> coordinate;
    std::optional<NastranBulkPosition> position;
};

struct NastranPositionLimits {
    structure::ExactLimits exact;
    std::uint64_t maximum_positions;
    std::uint64_t maximum_work;
    std::stop_token stop;
};

struct NastranPositionView {
    RootSnapshotLocation source;
    SourceDigest digest;
    SiLengthUnit length_unit;
    std::span<const NastranGridEnclosure> grids;
    std::uint64_t work;
};

// The endpoints are owned exact dyadics enclosing a position. No interval,
// including a singleton, issues the exact planar correspondence certificate.
class NastranPositionEnclosures {
public:
    NastranPositionEnclosures(NastranPositionEnclosures&&) noexcept;
    NastranPositionEnclosures& operator=(NastranPositionEnclosures&&) = delete;

    [[nodiscard]] static core::Result<NastranPositionEnclosures>
    prepare(const NastranCoordinateCards&, SiLengthUnit, std::span<const NastranGridId>,
            const NastranPositionLimits&, core::ResourcePlan&,
            std::optional<NastranPositionDiagnostic>* = nullptr) noexcept;

    [[nodiscard]] static core::Result<NastranPositionEnclosures>
    prepare(const NastranShellCards&, SiLengthUnit, std::span<const NastranGridId>,
            const NastranPositionLimits&, core::ResourcePlan&,
            std::optional<NastranPositionDiagnostic>* = nullptr) noexcept;

    [[nodiscard]] static core::Result<NastranPositionEnclosures>
    prepare(const NastranIsotropicShellCards&, SiLengthUnit, std::span<const NastranGridId>,
            const NastranPositionLimits&, core::ResourcePlan&,
            std::optional<NastranPositionDiagnostic>* = nullptr) noexcept;

    [[nodiscard]] std::optional<NastranPositionView> view() const& noexcept;
    [[nodiscard]] std::optional<NastranPositionView> view() const&& = delete;

private:
    NastranPositionEnclosures(RootSnapshotLocation, SourceDigest, SiLengthUnit,
                              core::PlannedBuffer<NastranGridEnclosure>&&, std::uint64_t) noexcept;
    [[nodiscard]] static core::Result<NastranPositionEnclosures>
    prepare_view(const NastranCoordinateView*, SiLengthUnit, std::span<const NastranGridId>,
                 const NastranPositionLimits&, core::ResourcePlan&,
                 std::optional<NastranPositionDiagnostic>*) noexcept;
    std::optional<RootSnapshotLocation> source_;
    SourceDigest digest_;
    SiLengthUnit length_unit_;
    core::PlannedBuffer<NastranGridEnclosure> grids_;
    std::uint64_t work_;
};

struct NastranSpcDirectionLimits {
    structure::ExactLimits exact;
    std::uint64_t maximum_work = 0;
    std::stop_token stop;
};

struct NastranSpcDirectionOrigin {
    NastranSupportOrigin source;
    NastranSupportSetId set{0};
    NastranGridId grid{0};
    NastranCoordinateId cd{0};
    NastranBulkPosition grid_field, component_field, card;
    std::optional<NastranBulkPosition> cd_field, motion_field;
    std::optional<NastranRealField> prescribed_motion;
};

enum class NastranSpcComponentKind : std::uint8_t { Translation, Rotation };

using NastranComponentAxis = std::array<NastranCartesianInterval, 3>;

struct NastranGridDisplacementAxes {
    std::array<NastranComponentAxis, 3> directions;
    std::uint64_t work = 0;
};

// The three low bits identify the GRID axes needed by a vector. Only those
// directions are readable by the caller. A sole cylindrical axial request
// remains defined on the polar axis.
[[nodiscard]] core::Result<NastranGridDisplacementAxes>
prepare_nastran_displacement_axes(const NastranIsotropicShellCards&, NastranGridId,
                                  std::uint8_t required_axes, const NastranSpcDirectionLimits&,
                                  core::ResourcePlan&,
                                  std::optional<NastranPositionDiagnostic>* = nullptr) noexcept;

struct NastranSpcComponentDirection {
    NastranSpcComponentKind kind;
    const NastranComponentAxis& basic;
};

// One selected source SPC/SPC1 term and GRID, including an existing THRU
// member. Returns unit directions in basic coordinates for selected DOFs.
// Term selection does not perform Case Control or assemble constraints.
class NastranSpcDirections {
public:
    NastranSpcDirections(NastranSpcDirections&&) noexcept;
    NastranSpcDirections& operator=(NastranSpcDirections&&) = delete;

    [[nodiscard]] static core::Result<NastranSpcDirections>
    prepare(const NastranCoordinateCards&, const NastranSupportCards&, std::size_t term,
            NastranGridId, const NastranSpcDirectionLimits&, core::ResourcePlan&,
            std::optional<NastranPositionDiagnostic>* = nullptr) noexcept;

    [[nodiscard]] const NastranSpcDirectionOrigin* origin() const& noexcept;
    [[nodiscard]] const NastranSpcDirectionOrigin* origin() const&& = delete;
    [[nodiscard]] std::optional<NastranSpcComponentDirection>
    direction(unsigned component) const& noexcept;
    [[nodiscard]] std::optional<NastranSpcComponentDirection> direction(unsigned) const&& = delete;
    [[nodiscard]] std::uint64_t work() const noexcept;

private:
    NastranSpcDirections(const NastranSpcDirectionOrigin&, std::uint8_t,
                         std::array<NastranComponentAxis, 3>&&, std::uint64_t) noexcept;
    std::optional<NastranSpcDirectionOrigin> origin_;
    std::uint8_t components_;
    std::array<NastranComponentAxis, 3> axes_;
    std::uint64_t work_;
};

} // namespace uvlas::io::detail
