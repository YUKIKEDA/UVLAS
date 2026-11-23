#include "checkpoint_manifest.hpp"

#include "json_object_reader.hpp"

#include <algorithm>
#include <charconv>
#include <limits>
#include <utility>

namespace uvlas::io::detail {
namespace {

constexpr auto recover = core::Operation::Recover;
constexpr auto save = core::Operation::SaveCheckpoint;
constexpr auto corrupt = core::FailureCode::CheckpointCorrupt;
constexpr auto incompatible = core::FailureCode::CheckpointIncompatible;
constexpr auto capacity = core::FailureCode::ResourceExceeded;
constexpr auto cancelled = core::FailureCode::Cancelled;
constexpr auto invalid = core::FailureCode::InvalidInput;
constexpr std::size_t name_capacity = 64;
constexpr std::size_t word_bytes = 8;
constexpr std::string_view schema = "uvlas.checkpoint/1";
constexpr std::string_view byte_order = "little";
constexpr std::string_view hex_digits = "0123456789abcdef";

struct FingerprintField {
    std::string_view name;
    SourceDigest CheckpointFingerprints::* digest;
};
constexpr std::array fingerprints{
    FingerprintField{"model", &CheckpointFingerprints::model},
    FingerprintField{"basis", &CheckpointFingerprints::basis},
    FingerprintField{"material", &CheckpointFingerprints::material},
    FingerprintField{"constraint", &CheckpointFingerprints::constraint},
    FingerprintField{"algorithm", &CheckpointFingerprints::algorithm},
    FingerprintField{"precision", &CheckpointFingerprints::precision},
    FingerprintField{"version", &CheckpointFingerprints::version},
};

struct ChunkEntry {
    std::array<char, name_capacity> name{};
    std::size_t name_size = 0;
    CheckpointDtype dtype = CheckpointDtype::F64;
    core::PlannedBuffer<core::Index> shape;
    std::uint64_t bytes = 0;
    SourceDigest digest{};

    [[nodiscard]] std::string_view name_view() const noexcept {
        return {name.data(), name_size};
    }
};

struct Scope {
    core::Operation operation;
    core::FailureCode malformed;
    const std::stop_token& stop;

    [[noreturn]] void fail(core::FailureCode code, core::FailureDetails details = {}) const {
        throw core::Failure{code, operation, details};
    }

    void require(bool condition, core::FailureCode code) const {
        if (!condition)
            fail(code);
    }

    void poll() const {
        require(!stop.stop_requested(), cancelled);
    }

    template <class T>
    [[nodiscard]] T get(core::Result<T> result) const {
        if (!result) {
            const auto& error = result.error();
            fail(error.code == invalid ? malformed : error.code, error.details);
        }
        return std::move(*result);
    }

    template <class T>
    [[nodiscard]] T get(std::optional<T> result) const {
        if (!result)
            fail(malformed);
        return std::move(*result);
    }
};

template <class T, class Call>
[[nodiscard]] core::Result<T> guarded(const Scope& scope, Call call) noexcept {
    return core::guard_entry(scope.operation, scope.malformed, {}, [&]() -> core::Result<T> {
        try {
            scope.poll();
            return call();
        } catch (const core::Failure& failure) {
            return std::unexpected(failure);
        }
    });
}

[[nodiscard]] bool lowercase_letter(char byte) noexcept {
    return byte >= 'a' && byte <= 'z';
}

[[nodiscard]] bool valid_name(std::string_view name) noexcept {
    return !name.empty() && name.size() <= name_capacity && lowercase_letter(name.front()) &&
           std::ranges::all_of(name, [](char byte) {
               return lowercase_letter(byte) || (byte >= '0' && byte <= '9') || byte == '_' ||
                      byte == '-';
           });
}

[[nodiscard]] std::string_view dtype_name(CheckpointDtype dtype) noexcept {
    switch (dtype) {
    case CheckpointDtype::F64:
        return "f64";
    case CheckpointDtype::I64:
        return "i64";
    case CheckpointDtype::U64:
        return "u64";
    }
    return {};
}

[[nodiscard]] CheckpointDtype read_dtype(std::string_view name, const Scope& scope) {
    for (const auto dtype : {CheckpointDtype::F64, CheckpointDtype::I64, CheckpointDtype::U64}) {
        if (name == dtype_name(dtype))
            return dtype;
    }
    scope.fail(incompatible);
}

[[nodiscard]] SourceDigest read_digest(std::string_view text, const Scope& scope) {
    SourceDigest result{};
    scope.require(text.size() == result.bytes.size() * 2, corrupt);
    for (std::size_t i = 0; i < result.bytes.size(); ++i) {
        const auto high = hex_digits.find(text[i * 2]);
        const auto low = hex_digits.find(text[i * 2 + 1]);
        scope.require(high != hex_digits.npos && low != hex_digits.npos, corrupt);
        result.bytes[i] = static_cast<std::byte>((high << 4) | low);
    }
    return result;
}

[[nodiscard]] std::string_view text_field(const InputObjectView& object, std::string_view key,
                                          const Scope& scope) {
    const auto node = scope.get(object.find(key));
    return scope.get(node.string());
}

[[nodiscard]] std::uint64_t unsigned_number(const InputNodeView& node, const Scope& scope) {
    if (const auto value = node.unsigned_integer())
        return *value;
    if (const auto value = node.signed_integer(); value && *value >= 0)
        return static_cast<std::uint64_t>(*value);
    scope.fail(corrupt);
}

void check_shape(const CheckpointChunkHeader& header, CheckpointChunkLimits limits,
                 const Scope& scope) {
    scope.require(header.shape.size() <= limits.maximum_rank, capacity);
    scope.require(header.byte_count <= limits.maximum_payload_bytes, capacity);
    const auto count = scope.get(core::checked_extent_product(header.shape, scope.operation));
    const auto bytes = scope.get(core::checked_byte_count(count, word_bytes, scope.operation));
    scope.require(std::in_range<core::Index>(bytes), capacity);
    scope.require(bytes == header.byte_count, scope.malformed);
}

void account_bytes(std::uint64_t add, std::uint64_t& total, std::uint64_t limit,
                   const Scope& scope) {
    scope.require(add <= limit - total, capacity);
    total += add;
}

[[nodiscard]] core::Index count_chunks(std::uint64_t count, const CheckpointManifestLimits& limits,
                                       const Scope& scope) {
    scope.require(count != 0, scope.malformed);
    scope.require(count <= limits.maximum_chunks && std::in_range<core::Index>(count), capacity);
    return static_cast<core::Index>(count);
}

template <class T, class Name>
void sort_and_check_names(std::span<T> rows, Name name, const Scope& scope) {
    std::ranges::sort(rows, {}, name);
    for (std::size_t i = 1; i < rows.size(); ++i) {
        scope.poll();
        scope.require(std::invoke(name, rows[i - 1]) != std::invoke(name, rows[i]),
                      scope.malformed);
    }
}

[[nodiscard]] ChunkEntry read_chunk(const InputNodeView& node,
                                    const CheckpointManifestLimits& limits,
                                    core::ResourcePlan& plan, const Scope& scope) {
    const auto object = scope.get(node.object());
    scope.require(object.size() == 5, corrupt);
    const auto name = text_field(object, "name", scope);
    scope.require(valid_name(name), corrupt);
    ChunkEntry result;
    std::ranges::copy(name, result.name.begin());
    result.name_size = name.size();
    result.dtype = read_dtype(text_field(object, "dtype", scope), scope);
    result.bytes = unsigned_number(scope.get(object.find("bytes")), scope);
    result.digest = read_digest(text_field(object, "sha256", scope), scope);
    const auto shape_node = scope.get(object.find("shape"));
    const auto shape = scope.get(shape_node.array());
    scope.require(static_cast<std::uint64_t>(shape.size()) <= limits.chunk.maximum_rank, capacity);
    result.shape = scope.get(plan.allocate<core::Index>(shape.size(), recover));
    for (core::Index i = 0; i < shape.size(); ++i) {
        scope.poll();
        const auto extent = unsigned_number(scope.get(shape.at(i)), scope);
        scope.require(std::in_range<core::Index>(extent), capacity);
        result.shape.values()[static_cast<std::size_t>(i)] = static_cast<core::Index>(extent);
    }
    check_shape({result.dtype, result.shape.values(), result.bytes, result.digest}, limits.chunk,
                scope);
    return result;
}

struct WriteRow {
    std::string_view name;
    CheckpointChunkHeader header{};
};

class TextOutput {
public:
    TextOutput(std::uint64_t limit, const Scope& scope, std::span<char> buffer = {}) noexcept
        : limit_(limit), scope_(scope), buffer_(buffer) {}

    void append(std::string_view text) {
        const auto start = size_;
        account_bytes(text.size(), size_, limit_, scope_);
        if (!buffer_.empty())
            std::ranges::copy(text, buffer_.subspan(static_cast<std::size_t>(start)).begin());
    }

    void quoted(std::string_view text) {
        append("\"");
        append(text);
        append("\"");
    }

    void number(std::uint64_t value) {
        std::array<char, std::numeric_limits<std::uint64_t>::digits10 + 1> storage{};
        const auto [end, error] =
            std::to_chars(storage.data(), storage.data() + storage.size(), value);
        scope_.require(error == std::errc{}, capacity);
        append({storage.data(), static_cast<std::size_t>(end - storage.data())});
    }

    void digest(const SourceDigest& value) {
        std::array<char, 64> text{};
        for (std::size_t i = 0; i < value.bytes.size(); ++i) {
            const auto byte = std::to_integer<unsigned>(value.bytes[i]);
            text[i * 2] = hex_digits[byte >> 4];
            text[i * 2 + 1] = hex_digits[byte & 0x0fU];
        }
        quoted({text.data(), text.size()});
    }

    [[nodiscard]] std::uint64_t size() const noexcept {
        return size_;
    }

private:
    std::uint64_t limit_;
    const Scope& scope_;
    std::span<char> buffer_;
    std::uint64_t size_ = 0;
};

void write_json(TextOutput& output, std::span<const WriteRow> rows,
                const CheckpointFingerprints& values, const Scope& scope) {
    output.append("{\"schema\":");
    output.quoted(schema);
    output.append(",\"byte_order\":");
    output.quoted(byte_order);
    output.append(",\"fingerprints\":{");
    for (std::size_t i = 0; i < fingerprints.size(); ++i) {
        if (i != 0)
            output.append(",");
        const auto& field = fingerprints[i];
        output.quoted(field.name);
        output.append(":");
        output.digest(values.*field.digest);
    }
    output.append("},\"chunks\":[");
    for (std::size_t i = 0; i < rows.size(); ++i) {
        scope.poll();
        const auto& row = rows[i];
        if (i != 0)
            output.append(",");
        output.append("{\"name\":");
        output.quoted(row.name);
        output.append(",\"dtype\":");
        output.quoted(dtype_name(row.header.dtype));
        output.append(",\"shape\":[");
        for (std::size_t d = 0; d < row.header.shape.size(); ++d) {
            scope.poll();
            if (d != 0)
                output.append(",");
            output.number(static_cast<std::uint64_t>(row.header.shape[d]));
        }
        output.append("],\"bytes\":");
        output.number(row.header.byte_count);
        output.append(",\"sha256\":");
        output.digest(row.header.digest);
        output.append("}");
    }
    output.append("]}\n");
}

} // namespace

struct CheckpointManifestStorage {
    core::PlannedBuffer<ChunkEntry> chunks;
};

CheckpointManifest::CheckpointManifest(
    core::PlannedBuffer<CheckpointManifestStorage> storage) noexcept
    : storage_(std::move(storage)) {}
CheckpointManifest::CheckpointManifest(CheckpointManifest&&) noexcept = default;
CheckpointManifest& CheckpointManifest::operator=(CheckpointManifest&&) noexcept = default;
CheckpointManifest::~CheckpointManifest() = default;

core::Index CheckpointManifest::size() const noexcept {
    return storage_.size() == 0 ? 0 : static_cast<core::Index>(storage_.values()[0].chunks.size());
}

core::Result<CheckpointManifestEntry> CheckpointManifest::at(core::Index index) const& noexcept {
    if (index < 0 || index >= size())
        return std::unexpected(core::Failure{invalid, recover, {}});
    const auto& chunk = storage_.values()[0].chunks.values()[static_cast<std::size_t>(index)];
    return CheckpointManifestEntry{chunk.name_view(),
                                   {chunk.dtype, chunk.shape.values(), chunk.bytes, chunk.digest}};
}

core::Result<CheckpointManifest> CheckpointManifest::read(std::span<const std::byte> bytes,
                                                          SourceId source,
                                                          const CheckpointFingerprints& expected,
                                                          const CheckpointManifestLimits& limits,
                                                          core::ResourcePlan& plan,
                                                          const std::stop_token& stop) noexcept {
    const Scope scope{recover, corrupt, stop};
    return guarded<CheckpointManifest>(scope, [&] {
        std::optional<InputDiagnostic> diagnostic;
        auto document = read_json_object(bytes, source, limits.json, plan, diagnostic);
        if (!document) {
            const auto& error = document.error();
            scope.fail(error.code == invalid || error.code == core::FailureCode::Nonfinite
                           ? corrupt
                           : error.code,
                       error.details);
        }
        scope.poll();
        const auto view = scope.get(document->view());
        const auto root = scope.get(view.root.object());
        scope.require(text_field(root, "schema", scope) == schema, incompatible);
        scope.require(text_field(root, "byte_order", scope) == byte_order, incompatible);
        scope.require(root.size() == 4, corrupt);
        const auto fingerprint_node = scope.get(root.find("fingerprints"));
        const auto fingerprint_object = scope.get(fingerprint_node.object());
        scope.require(fingerprint_object.size() == static_cast<core::Index>(fingerprints.size()),
                      corrupt);
        CheckpointFingerprints actual{};
        for (const auto& field : fingerprints)
            actual.*field.digest =
                read_digest(text_field(fingerprint_object, field.name, scope), scope);
        scope.require(actual == expected, incompatible);
        const auto series = scope.get(root.find("chunks"));
        const auto chunks = scope.get(series.array());
        const auto count = count_chunks(static_cast<std::uint64_t>(chunks.size()), limits, scope);
        auto entries = scope.get(plan.allocate<ChunkEntry>(count, recover));
        std::uint64_t payload = 0;
        for (core::Index i = 0; i < count; ++i) {
            scope.poll();
            auto entry = read_chunk(scope.get(chunks.at(i)), limits, plan, scope);
            account_bytes(entry.bytes, payload, limits.maximum_total_payload_bytes, scope);
            entries.values()[static_cast<std::size_t>(i)] = std::move(entry);
        }
        sort_and_check_names(entries.values(), &ChunkEntry::name_view, scope);
        auto storage = scope.get(plan.allocate<CheckpointManifestStorage>(1, recover));
        storage.values()[0].chunks = std::move(entries);
        scope.poll();
        return CheckpointManifest{std::move(storage)};
    });
}

core::Result<core::PlannedBuffer<char>>
write_checkpoint_manifest(std::span<const NamedCheckpointChunk> chunks,
                          const CheckpointFingerprints& expected,
                          const CheckpointManifestLimits& limits, core::ResourcePlan& plan,
                          const std::stop_token& stop) noexcept {
    const Scope scope{save, invalid, stop};
    return guarded<core::PlannedBuffer<char>>(scope, [&] {
        const auto count = count_chunks(chunks.size(), limits, scope);
        auto rows = scope.get(plan.allocate<WriteRow>(count, save));
        std::uint64_t payload = 0;
        for (std::size_t i = 0; i < chunks.size(); ++i) {
            scope.poll();
            const auto& entry = chunks[i];
            scope.require(valid_name(entry.name), invalid);
            const auto view = scope.get(entry.chunk.get().view());
            check_shape(view.header, limits.chunk, scope);
            account_bytes(view.header.byte_count, payload, limits.maximum_total_payload_bytes,
                          scope);
            rows.values()[i] = {entry.name, view.header};
        }
        sort_and_check_names(rows.values(), &WriteRow::name, scope);
        TextOutput measure{limits.json.maximum_source_bytes, scope};
        write_json(measure, rows.values(), expected, scope);
        scope.require(std::in_range<core::Index>(measure.size()), capacity);
        auto buffer =
            scope.get(plan.allocate<char>(static_cast<core::Index>(measure.size()), save));
        TextOutput write{measure.size(), scope, buffer.values()};
        write_json(write, rows.values(), expected, scope);
        scope.poll();
        return buffer;
    });
}

CheckpointPayloadSet::CheckpointPayloadSet(
    core::PlannedBuffer<std::optional<CheckpointChunk>> chunks) noexcept
    : chunks_(std::move(chunks)) {}

core::Index CheckpointPayloadSet::size() const noexcept {
    return static_cast<core::Index>(chunks_.size());
}

core::Result<std::reference_wrapper<const CheckpointChunk>>
CheckpointPayloadSet::at(core::Index index) const& noexcept {
    if (index < 0 || index >= size())
        return std::unexpected(core::Failure{invalid, recover, {}});
    const auto& chunk = chunks_.values()[static_cast<std::size_t>(index)];
    if (!chunk)
        return std::unexpected(core::Failure{corrupt, recover, {}});
    return std::cref(*chunk);
}

core::Result<CheckpointPayloadSet> authenticate_checkpoint_payloads(
    const CheckpointManifest& manifest, std::span<const CheckpointPayload> payloads,
    CheckpointChunkLimits limits, core::ResourcePlan& plan, const std::stop_token& stop) noexcept {
    const Scope scope{recover, corrupt, stop};
    return guarded<CheckpointPayloadSet>(scope, [&] {
        const auto count = manifest.size();
        scope.require(count > 0 && static_cast<std::uint64_t>(count) == payloads.size(), corrupt);
        auto sorted = scope.get(plan.allocate<CheckpointPayload>(count, recover));
        for (std::size_t i = 0; i < payloads.size(); ++i) {
            scope.poll();
            scope.require(valid_name(payloads[i].name), corrupt);
            sorted.values()[i] = payloads[i];
        }
        sort_and_check_names(sorted.values(), &CheckpointPayload::name, scope);
        auto result = scope.get(plan.allocate<std::optional<CheckpointChunk>>(count, recover));
        for (core::Index i = 0; i < count; ++i) {
            scope.poll();
            const auto entry = scope.get(manifest.at(i));
            const auto index = static_cast<std::size_t>(i);
            const auto& observed = sorted.values()[index];
            scope.require(observed.name == entry.name, corrupt);
            result.values()[index].emplace(scope.get(
                CheckpointChunk::authenticate({entry.header, observed.bytes}, limits, plan, stop)));
        }
        scope.poll();
        return CheckpointPayloadSet{std::move(result)};
    });
}

} // namespace uvlas::io::detail
