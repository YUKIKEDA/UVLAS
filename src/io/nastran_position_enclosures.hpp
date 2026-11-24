#pragma once

#include "nastran_shell_cards.hpp"

#include <uvlas/io/si_case.hpp>

namespace uvlas::io::detail {

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

} // namespace uvlas::io::detail
