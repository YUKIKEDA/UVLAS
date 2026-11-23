#pragma once

#include <optional>
#include <stop_token>
#include <uvlas/core/resource_plan.hpp>
#include <uvlas/io/input_document.hpp>
#include <variant>

namespace uvlas::io::detail {

class ReferenceFileSnapshot;

enum class CheckpointDtype : std::uint8_t { F64, I64, U64 };

using CheckpointInput = std::variant<std::span<const double>, std::span<const std::int64_t>,
                                     std::span<const std::uint64_t>>;
using CheckpointValues =
    std::variant<core::PlannedBuffer<double>, core::PlannedBuffer<std::int64_t>,
                 core::PlannedBuffer<std::uint64_t>>;

struct CheckpointChunkHeader {
    CheckpointDtype dtype = CheckpointDtype::F64;
    std::span<const core::Index> shape;
    std::uint64_t byte_count = 0;
    SourceDigest digest{};
};

struct CheckpointChunkView {
    CheckpointChunkHeader header;
    std::span<const std::byte> bytes;
};

struct CheckpointChunkLimits {
    std::uint64_t maximum_rank;
    std::uint64_t maximum_payload_bytes;
};

/// Owns one authenticated chunk. The caller supplies the manifest fields and keeps
/// borrowed inputs stable for the duration of each call.
class CheckpointChunk {
public:
    CheckpointChunk(CheckpointChunk&& other) noexcept;
    CheckpointChunk& operator=(CheckpointChunk&& other) noexcept;

    [[nodiscard]] static core::Result<CheckpointChunk>
    encode(std::span<const core::Index> shape, CheckpointInput values, CheckpointChunkLimits limits,
           core::ResourcePlan& plan, const std::stop_token& stop) noexcept;

    [[nodiscard]] static core::Result<CheckpointChunk>
    authenticate(CheckpointChunkView input, CheckpointChunkLimits limits, core::ResourcePlan& plan,
                 const std::stop_token& stop) noexcept;

    /// Reuse bytes already read and hashed from a stable native file handle. The
    /// snapshot is consumed only after shape, size and digest agree.
    [[nodiscard]] static core::Result<CheckpointChunk>
    from_verified_file(CheckpointChunkHeader, ReferenceFileSnapshot&&, CheckpointChunkLimits,
                       core::ResourcePlan&, const std::stop_token&) noexcept;

    /// Borrow until destruction or ownership transfer. A moved-from chunk has no view.
    [[nodiscard]] std::optional<CheckpointChunkView> view() const& noexcept;
    [[nodiscard]] std::optional<CheckpointChunkView> view() const&& = delete;

    /// Materialize an independent typed array with the original 64-bit words.
    [[nodiscard]] core::Result<CheckpointValues> decode(core::ResourcePlan& plan,
                                                        const std::stop_token& stop) const noexcept;

private:
    template <class T>
    [[nodiscard]] static core::Result<CheckpointChunk>
    encode_typed(std::span<const core::Index> shape, std::span<const T> values,
                 CheckpointDtype dtype, CheckpointChunkLimits limits, core::ResourcePlan& plan,
                 const std::stop_token& stop) noexcept;

    CheckpointChunk(CheckpointDtype dtype, SourceDigest digest, core::Index element_count,
                    core::PlannedBuffer<core::Index> shape,
                    core::PlannedBuffer<std::byte> bytes) noexcept;

    std::optional<CheckpointDtype> dtype_;
    SourceDigest digest_;
    core::Index element_count_;
    core::PlannedBuffer<core::Index> shape_;
    core::PlannedBuffer<std::byte> bytes_;
};

} // namespace uvlas::io::detail
