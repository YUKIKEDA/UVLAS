#pragma once

#include "aero/pchip_rational_slots.hpp"
#include "case_motion_rotation_private.hpp"
#include "io/si_exact_words.hpp"

#include <uvlas/io/si_case.hpp>

namespace uvlas::coupling::detail {

using CaseMotionField = decltype(io::schema::AeroAttachmentPrescribed::motion);

/// Fixed stack-payload envelope for the interval imports. Exact rational heap
/// admission is checked separately by the adopted PG and SI arithmetic.
inline constexpr std::uint64_t case_motion_frame_interval_scratch_bytes =
    12 * sizeof(aero::detail::ExactInteger) + 48 * sizeof(motion::m::Integer);

struct CaseMotionFrameLimits {
    structure::ExactLimits exact;
    std::uint64_t maximum_work{}; ///< Recipe inspections, algebra and word/interval work.
    std::uint64_t
        maximum_scratch_bytes{}; ///< Fixed interval-import payload, checked at preparation.
};

struct CaseMotionFrameView {
    const motion::Vector& origin;   ///< World coordinates, in metres.
    const motion::Matrix& rotation; ///< Row-major Q0. Columns are world images of the local axes.
};

struct CaseMotionFrameInfo {
    io::SchemaIdentity source;
    std::uint64_t origin_component_node; ///< First original numeric slot identifies the Motion.
    std::uint64_t work;
};

/// Borrowed save of the three world-origin and nine row-major rotation components.
struct CaseMotionFrameSnapshotView {
    CaseMotionFrameSnapshotView(io::SchemaIdentity identity, std::uint64_t origin_node,
                                std::array<aero::detail::pchip::RationalView, 12> values) noexcept
        : source(identity), origin_component_node(origin_node), components(values) {}

    io::SchemaIdentity source;
    std::uint64_t origin_component_node{};
    std::array<aero::detail::pchip::RationalView, 12> components;
};

class CaseMotionFrameSnapshot {
public:
    CaseMotionFrameSnapshot(CaseMotionFrameSnapshot&&) noexcept = default;
    CaseMotionFrameSnapshot& operator=(CaseMotionFrameSnapshot&&) noexcept = default;
    CaseMotionFrameSnapshot(const CaseMotionFrameSnapshot&) = delete;
    CaseMotionFrameSnapshot& operator=(const CaseMotionFrameSnapshot&) = delete;

    [[nodiscard]] core::Result<CaseMotionFrameSnapshotView> view() const& noexcept;
    [[nodiscard]] core::Result<CaseMotionFrameSnapshotView> view() const&& = delete;

private:
    friend class CaseMotionFrame;
    CaseMotionFrameSnapshot(core::PlannedBuffer<io::detail::SiRationalWords>,
                            CaseMotionFrameInfo) noexcept;
    core::PlannedBuffer<io::detail::SiRationalWords> components_;
    CaseMotionFrameInfo info_;
};

/// Owns all twelve exact SI values, their interval enclosures and motion provenance.
class CaseMotionFrame {
public:
    CaseMotionFrame(CaseMotionFrame&&) noexcept;
    CaseMotionFrame& operator=(CaseMotionFrame&&) noexcept;
    CaseMotionFrame(const CaseMotionFrame&) = delete;
    CaseMotionFrame& operator=(const CaseMotionFrame&) = delete;
    ~CaseMotionFrame();

    /// The field must come from this SiCase's immutable typed source.
    [[nodiscard]] static core::Result<CaseMotionFrame>
    prepare(const io::SiCase&, CaseMotionField, const CaseMotionFrameLimits&, core::ResourcePlan&,
            std::stop_token = {}, core::FailureDetails = {}) noexcept;

    /// Match the original typed node and compare all twelve exact SI values before use.
    [[nodiscard]] core::Result<void> validate_for(const io::SiCase&, CaseMotionField,
                                                  const CaseMotionFrameLimits&, core::ResourcePlan&,
                                                  std::stop_token = {},
                                                  core::FailureDetails = {}) const noexcept;

    /// Compare separately owned, canonical values without allocation. Preparation bounds
    /// the integer widths; cancellation is checked before each rational comparison.
    [[nodiscard]] core::Result<void> match(const CaseMotionFrame&, std::stop_token = {},
                                           core::FailureDetails = {}) const noexcept;

    [[nodiscard]] core::Result<CaseMotionFrameView> view() const& noexcept;
    [[nodiscard]] core::Result<CaseMotionFrameView> view() const&& = delete;
    [[nodiscard]] core::Result<CaseMotionFrameInfo> info() const noexcept;
    [[nodiscard]] core::Result<CaseMotionFrameSnapshot>
    snapshot(const CaseMotionFrameLimits&, core::ResourcePlan&, const std::stop_token& = {},
             core::FailureDetails = {}) const noexcept;
    [[nodiscard]] core::Result<CaseMotionFrame>
    restore(const CaseMotionFrameSnapshotView&, const CaseMotionFrameLimits&, core::ResourcePlan&,
            const std::stop_token& = {}, core::FailureDetails = {}) const noexcept;

private:
    struct Data;
    explicit CaseMotionFrame(core::PlannedBuffer<Data>, CaseMotionFrameInfo) noexcept;
    core::PlannedBuffer<Data> data_;
    CaseMotionFrameInfo info_;
};

} // namespace uvlas::coupling::detail
