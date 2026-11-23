#include "checkpoint_chunk.hpp"

#include "source_digest.hpp"

#include <algorithm>
#include <bit>
#include <limits>
#include <utility>

namespace uvlas::io::detail {
namespace {

constexpr auto save = core::Operation::SaveCheckpoint;
constexpr auto recover = core::Operation::Recover;
constexpr std::size_t word_bytes = 8;
constexpr auto block_bytes = std::size_t{64} * 1024;
constexpr core::Index block_words = block_bytes / word_bytes;
static_assert(sizeof(double) == word_bytes && std::numeric_limits<double>::is_iec559);
static_assert(sizeof(std::int64_t) == word_bytes && sizeof(std::uint64_t) == word_bytes);

struct ChunkStorage {
    core::Index words;
    core::Index bytes;
    core::Index rank;
};

[[nodiscard]] std::unexpected<core::Failure> fail(core::FailureCode code,
                                                  core::Operation operation) noexcept {
    return std::unexpected(core::Failure{code, operation, {}});
}

[[nodiscard]] core::Result<ChunkStorage> measure(std::span<const core::Index> shape,
                                                 CheckpointChunkLimits limits,
                                                 core::Operation operation,
                                                 const std::stop_token& stop) noexcept {
    if (stop.stop_requested()) {
        return fail(core::FailureCode::Cancelled, operation);
    }
    if (shape.size() > limits.maximum_rank || !std::in_range<core::Index>(shape.size())) {
        return fail(core::FailureCode::ResourceExceeded, operation);
    }
    const auto words = core::checked_extent_product(shape, operation);
    if (!words) {
        if (operation == recover && words.error().code == core::FailureCode::InvalidInput) {
            return fail(core::FailureCode::CheckpointCorrupt, operation);
        }
        return std::unexpected(words.error());
    }
    const auto bytes = core::checked_byte_count(*words, word_bytes, operation);
    if (!bytes) {
        return std::unexpected(bytes.error());
    }
    if (*bytes > limits.maximum_payload_bytes || !std::in_range<core::Index>(*bytes)) {
        return fail(core::FailureCode::ResourceExceeded, operation);
    }
    return ChunkStorage{*words, static_cast<core::Index>(*bytes),
                        static_cast<core::Index>(shape.size())};
}

[[nodiscard]] core::Result<void> admit_storage(const ChunkStorage& size,
                                               const core::ResourcePlan& plan,
                                               core::Operation operation) noexcept {
    const auto shape_bytes = core::checked_byte_count(size.rank, sizeof(core::Index), operation);
    if (!shape_bytes) {
        return std::unexpected(shape_bytes.error());
    }
    const auto payload = static_cast<std::uint64_t>(size.bytes);
    if (payload > plan.remaining_bytes() || *shape_bytes > plan.remaining_bytes() - payload) {
        return fail(core::FailureCode::ResourceExceeded, operation);
    }
    return {};
}

[[nodiscard]] core::Result<core::PlannedBuffer<core::Index>>
copy_shape(std::span<const core::Index> shape, core::Index rank, core::ResourcePlan& plan,
           core::Operation operation, const std::stop_token& stop) noexcept {
    if (stop.stop_requested()) {
        return fail(core::FailureCode::Cancelled, operation);
    }
    auto output = plan.allocate<core::Index>(rank, operation);
    if (!output) {
        return std::unexpected(output.error());
    }
    for (std::size_t offset = 0; offset < shape.size();) {
        if (stop.stop_requested()) {
            return fail(core::FailureCode::Cancelled, operation);
        }
        const auto batch = shape.subspan(
            offset, std::min(block_bytes / sizeof(core::Index), shape.size() - offset));
        std::ranges::copy(batch, output->values().subspan(offset).begin());
        offset += batch.size();
    }
    return output;
}

template <class T>
[[nodiscard]] core::Result<SourceDigest> write_words(std::span<const T> input,
                                                     std::span<std::byte> output,
                                                     const std::stop_token& stop) noexcept {
    boost::hash2::sha2_256 hash;
    while (!input.empty()) {
        if (stop.stop_requested()) {
            return fail(core::FailureCode::Cancelled, save);
        }
        const auto batch =
            input.first(std::min(input.size(), static_cast<std::size_t>(block_words)));
        for (std::size_t i = 0; i < batch.size(); ++i) {
            const auto word = std::bit_cast<std::uint64_t>(batch[i]);
            for (std::size_t b = 0; b < word_bytes; ++b) {
                output[i * word_bytes + b] = static_cast<std::byte>((word >> (8 * b)) & 0xffU);
            }
        }
        const auto encoded_bytes = batch.size() * word_bytes;
        hash.update(output.data(), encoded_bytes);
        input = input.subspan(batch.size());
        output = output.subspan(encoded_bytes);
    }
    if (stop.stop_requested()) {
        return fail(core::FailureCode::Cancelled, save);
    }
    return finish_source_digest(hash);
}

template <class T>
[[nodiscard]] core::Result<CheckpointValues> read_words(std::span<const std::byte> input,
                                                        core::Index count, core::ResourcePlan& plan,
                                                        const std::stop_token& stop) noexcept {
    auto values = plan.allocate<T>(count, recover);
    if (!values) {
        return std::unexpected(values.error());
    }
    auto output = values->values();
    while (!output.empty()) {
        if (stop.stop_requested()) {
            return fail(core::FailureCode::Cancelled, recover);
        }
        const auto batch =
            output.first(std::min(output.size(), static_cast<std::size_t>(block_words)));
        for (std::size_t i = 0; i < batch.size(); ++i) {
            std::uint64_t word = 0;
            for (std::size_t b = 0; b < word_bytes; ++b) {
                word |= std::to_integer<std::uint64_t>(input[i * word_bytes + b]) << (8 * b);
            }
            batch[i] = std::bit_cast<T>(word);
        }
        input = input.subspan(batch.size() * word_bytes);
        output = output.subspan(batch.size());
    }
    if (stop.stop_requested()) {
        return fail(core::FailureCode::Cancelled, recover);
    }
    return CheckpointValues{std::move(*values)};
}

} // namespace

CheckpointChunk::CheckpointChunk(CheckpointDtype dtype, SourceDigest digest,
                                 core::Index element_count, core::PlannedBuffer<core::Index> shape,
                                 core::PlannedBuffer<std::byte> bytes) noexcept
    : dtype_(dtype),
      digest_(digest),
      element_count_(element_count),
      shape_(std::move(shape)),
      bytes_(std::move(bytes)) {}

CheckpointChunk::CheckpointChunk(CheckpointChunk&& other) noexcept
    : dtype_(std::exchange(other.dtype_, std::nullopt)),
      digest_(other.digest_),
      element_count_(other.element_count_),
      shape_(std::move(other.shape_)),
      bytes_(std::move(other.bytes_)) {}

CheckpointChunk& CheckpointChunk::operator=(CheckpointChunk&& other) noexcept {
    if (this != &other) {
        dtype_ = std::exchange(other.dtype_, std::nullopt);
        digest_ = other.digest_;
        element_count_ = other.element_count_;
        shape_ = std::move(other.shape_);
        bytes_ = std::move(other.bytes_);
    }
    return *this;
}

core::Result<CheckpointChunk> CheckpointChunk::encode(std::span<const core::Index> shape,
                                                      CheckpointInput values,
                                                      CheckpointChunkLimits limits,
                                                      core::ResourcePlan& plan,
                                                      const std::stop_token& stop) noexcept {
    if (const auto* f64 = std::get_if<std::span<const double>>(&values)) {
        return encode_typed(shape, *f64, CheckpointDtype::F64, limits, plan, stop);
    }
    if (const auto* i64 = std::get_if<std::span<const std::int64_t>>(&values)) {
        return encode_typed(shape, *i64, CheckpointDtype::I64, limits, plan, stop);
    }
    if (const auto* u64 = std::get_if<std::span<const std::uint64_t>>(&values)) {
        return encode_typed(shape, *u64, CheckpointDtype::U64, limits, plan, stop);
    }
    return fail(core::FailureCode::InvalidInput, save);
}

template <class T>
core::Result<CheckpointChunk>
CheckpointChunk::encode_typed(std::span<const core::Index> shape, std::span<const T> values,
                              CheckpointDtype dtype, CheckpointChunkLimits limits,
                              core::ResourcePlan& plan, const std::stop_token& stop) noexcept {
    const auto size = measure(shape, limits, save, stop);
    if (!size) {
        return std::unexpected(size.error());
    }
    if (values.size() != static_cast<std::uint64_t>(size->words)) {
        return fail(core::FailureCode::InvalidInput, save);
    }
    const auto admitted = admit_storage(*size, plan, save);
    if (!admitted) {
        return std::unexpected(admitted.error());
    }
    auto owned_shape = copy_shape(shape, size->rank, plan, save, stop);
    if (!owned_shape) {
        return std::unexpected(owned_shape.error());
    }
    auto buffer = plan.allocate<std::byte>(size->bytes, save);
    if (!buffer) {
        return std::unexpected(buffer.error());
    }
    const auto digest = write_words(values, buffer->values(), stop);
    if (!digest) {
        return std::unexpected(digest.error());
    }
    return CheckpointChunk{dtype, *digest, size->words, std::move(*owned_shape),
                           std::move(*buffer)};
}

core::Result<CheckpointChunk> CheckpointChunk::authenticate(CheckpointChunkView input,
                                                            CheckpointChunkLimits limits,
                                                            core::ResourcePlan& plan,
                                                            const std::stop_token& stop) noexcept {
    const auto& header = input.header;
    if (header.dtype != CheckpointDtype::F64 && header.dtype != CheckpointDtype::I64 &&
        header.dtype != CheckpointDtype::U64) {
        return fail(core::FailureCode::CheckpointIncompatible, recover);
    }
    const auto size = measure(header.shape, limits, recover, stop);
    if (!size) {
        return std::unexpected(size.error());
    }
    if (static_cast<std::uint64_t>(size->bytes) != header.byte_count ||
        input.bytes.size() != header.byte_count) {
        return fail(core::FailureCode::CheckpointCorrupt, recover);
    }
    const auto admitted = admit_storage(*size, plan, recover);
    if (!admitted) {
        return std::unexpected(admitted.error());
    }
    boost::hash2::sha2_256 hash;
    for (std::size_t offset = 0; offset < input.bytes.size();) {
        if (stop.stop_requested()) {
            return fail(core::FailureCode::Cancelled, recover);
        }
        const auto batch =
            input.bytes.subspan(offset, std::min(block_bytes, input.bytes.size() - offset));
        hash.update(batch.data(), batch.size());
        offset += batch.size();
    }
    if (finish_source_digest(hash) != header.digest) {
        return fail(core::FailureCode::CheckpointCorrupt, recover);
    }
    auto owned_shape = copy_shape(header.shape, size->rank, plan, recover, stop);
    if (!owned_shape) {
        return std::unexpected(owned_shape.error());
    }
    auto bytes = plan.allocate<std::byte>(size->bytes, recover);
    if (!bytes) {
        return std::unexpected(bytes.error());
    }
    for (std::size_t offset = 0; offset < input.bytes.size();) {
        if (stop.stop_requested()) {
            return fail(core::FailureCode::Cancelled, recover);
        }
        const auto batch =
            input.bytes.subspan(offset, std::min(block_bytes, input.bytes.size() - offset));
        std::ranges::copy(batch, bytes->values().subspan(offset).begin());
        offset += batch.size();
    }
    if (stop.stop_requested()) {
        return fail(core::FailureCode::Cancelled, recover);
    }
    return CheckpointChunk{header.dtype, header.digest, size->words, std::move(*owned_shape),
                           std::move(*bytes)};
}

std::optional<CheckpointChunkView> CheckpointChunk::view() const& noexcept {
    if (!dtype_) {
        return std::nullopt;
    }
    return CheckpointChunkView{{*dtype_, shape_.values(), bytes_.size(), digest_}, bytes_.values()};
}

core::Result<CheckpointValues> CheckpointChunk::decode(core::ResourcePlan& plan,
                                                       const std::stop_token& stop) const noexcept {
    if (stop.stop_requested()) {
        return fail(core::FailureCode::Cancelled, recover);
    }
    if (!dtype_) {
        return fail(core::FailureCode::InvalidInput, recover);
    }
    switch (*dtype_) {
    case CheckpointDtype::F64:
        return read_words<double>(bytes_.values(), element_count_, plan, stop);
    case CheckpointDtype::I64:
        return read_words<std::int64_t>(bytes_.values(), element_count_, plan, stop);
    case CheckpointDtype::U64:
        return read_words<std::uint64_t>(bytes_.values(), element_count_, plan, stop);
    }
    return fail(core::FailureCode::CheckpointIncompatible, recover);
}

} // namespace uvlas::io::detail
