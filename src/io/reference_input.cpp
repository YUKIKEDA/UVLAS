#include "reference_input.hpp"

#include "case_schema_storage.hpp"
#include "utf8.hpp"

#include <algorithm>
#include <utility>

namespace uvlas::io::detail {
namespace {

constexpr auto operation = core::Operation::PrepareModel;

[[nodiscard]] auto fail(core::FailureCode code) noexcept {
    return std::unexpected(core::Failure{code, operation});
}

[[nodiscard]] auto invalid() noexcept {
    return fail(core::FailureCode::InvalidInput);
}

[[nodiscard]] bool cancelled(const std::stop_token& stop,
                             ReferenceInputDiagnostic* diagnostic) noexcept {
    if (!stop.stop_requested()) {
        return false;
    }
    if (diagnostic) {
        diagnostic->reset();
    }
    return true;
}

void mark(ReferenceInputDiagnostic* diagnostic, ReferenceInputStage stage) noexcept {
    if (diagnostic) {
        diagnostic->stage = stage;
    }
}

void clear(ReferenceInputDiagnostic* diagnostic) noexcept {
    if (diagnostic) {
        diagnostic->reset();
    }
}

[[nodiscard]] core::Result<InputDocument>
read_document(std::span<const std::byte> bytes, SourceId source, const InputLimits& limits,
              core::ResourcePlan& plan, const std::stop_token& stop,
              ReferenceInputDiagnostic* diagnostic) noexcept {
    mark(diagnostic, ReferenceInputStage::Document);
    std::optional<InputDiagnostic> local;
    auto document =
        read_input_document(bytes, source, limits, plan, diagnostic ? diagnostic->document : local);
    if (cancelled(stop, diagnostic)) {
        return fail(core::FailureCode::Cancelled);
    }
    return document;
}

[[nodiscard]] core::Result<SchemaCase>
read_case(const InputDocument& document, const SchemaLimits& limits, core::ResourcePlan& plan,
          const std::stop_token& stop, ReferenceInputDiagnostic* diagnostic) noexcept {
    mark(diagnostic, ReferenceInputStage::Schema);
    auto schema =
        validate_case_schema(document, limits, plan, diagnostic ? &diagnostic->schema : nullptr);
    if (cancelled(stop, diagnostic)) {
        return fail(core::FailureCode::Cancelled);
    }
    return schema;
}

#if defined(_WIN32) || defined(__linux__)
struct ManifestDocument {
    InputDocument document;
    ReferenceFileId identity;
};

[[nodiscard]] core::Result<ManifestDocument>
read_manifest_document(ManifestReader& reader, SourceId source, ReferenceFileLimits file_limits,
                       const InputLimits& document_limits, core::ResourcePlan& plan,
                       const std::stop_token& stop, ReferenceInputDiagnostic* diagnostic) noexcept {
    mark(diagnostic, ReferenceInputStage::ManifestRead);
    file_limits.maximum_file_bytes =
        std::min(file_limits.maximum_file_bytes, document_limits.maximum_source_bytes);
    const auto snapshot = read_manifest_file(reader, file_limits, plan, stop,
                                             diagnostic ? &diagnostic->file : nullptr);
    if (cancelled(stop, diagnostic)) {
        return fail(core::FailureCode::Cancelled);
    }
    if (!snapshot) {
        return std::unexpected(snapshot.error());
    }
    const auto verified = snapshot->view();
    if (!verified) {
        return invalid();
    }
    auto document = read_document(verified->bytes, source, document_limits, plan, stop, diagnostic);
    if (!document) {
        return std::unexpected(document.error());
    }
    mark(diagnostic, ReferenceInputStage::Identity);
    const auto parsed = document->view();
    if (!parsed || parsed->identity !=
                       SourceIdentity{source, verified->file_info.byte_count, verified->digest}) {
        return invalid();
    }
    // Destroy the temporary file payload before the schema makes its own source copy.
    return ManifestDocument{std::move(*document), verified->file_info.identity};
}
#endif

} // namespace

void ReferenceInputDiagnostic::reset() noexcept {
    stage.reset();
    path.reset();
    open.reset();
    file.reset();
    document.reset();
    schema.reset();
}

ReferenceInput::ReferenceInput(SchemaCase schema, Origin origin) noexcept
    : schema_(std::move(schema)), origin_(std::move(origin)) {}

core::Result<ReferenceInputView> ReferenceInput::view() const& noexcept {
    if (!schema_.identity()) {
        return invalid();
    }
    if (const auto* held = std::get_if<FileOrigin>(&origin_)) {
        return ReferenceInputView{schema_, &held->root, held->identity};
    }

    if (const auto* memory = std::get_if<MemoryOrigin>(&origin_)) {
        return ReferenceInputView{schema_, memory->root ? &*memory->root : nullptr, std::nullopt};
    }
    return invalid();
}

core::Result<ReferenceInput> ReferenceInput::from_file(
    const std::filesystem::path& parent, std::string_view filename, SourceId source,
    ManifestLimits manifest_limits, const ReferenceInputLimits& limits, core::ResourcePlan& plan,
    const std::stop_token& stop, ReferenceInputDiagnostic* diagnostic) noexcept {
    clear(diagnostic);
    if (stop.stop_requested()) {
        return fail(core::FailureCode::Cancelled);
    }
#if defined(_WIN32) || defined(__linux__)
    mark(diagnostic, ReferenceInputStage::ManifestName);
    const auto name_limit = manifest_limits.maximum_name_bytes;
    if (filename.size() <= name_limit &&
        (filename.find_first_of("/\\") != std::string_view::npos || !valid_utf8(filename))) {
        return invalid();
    }
    auto name = normalize_reference_path(filename, {name_limit, name_limit, 1}, plan, stop,
                                         diagnostic ? &diagnostic->path : nullptr);
    if (cancelled(stop, diagnostic)) {
        return fail(core::FailureCode::Cancelled);
    }
    if (!name) {
        return std::unexpected(name.error());
    }

    mark(diagnostic, ReferenceInputStage::Root);
    auto root = ReadRoot::open(parent, stop);
    if (cancelled(stop, diagnostic)) {
        return fail(core::FailureCode::Cancelled);
    }
    if (!root) {
        return std::unexpected(root.error());
    }

    mark(diagnostic, ReferenceInputStage::ManifestOpen);
    auto reader = root->open_file(*name, plan, stop, diagnostic ? &diagnostic->open : nullptr);
    if (cancelled(stop, diagnostic)) {
        return fail(core::FailureCode::Cancelled);
    }
    if (!reader) {
        return std::unexpected(reader.error());
    }

    const auto manifest = read_manifest_document(*reader, source, manifest_limits.file,
                                                 limits.document, plan, stop, diagnostic);
    if (!manifest) {
        return std::unexpected(manifest.error());
    }

    auto schema = read_case(manifest->document, limits.schema, plan, stop, diagnostic);
    if (!schema) {
        return std::unexpected(schema.error());
    }

    if (cancelled(stop, diagnostic)) {
        return fail(core::FailureCode::Cancelled);
    }

    clear(diagnostic);
    return ReferenceInput{std::move(*schema),
                          FileOrigin{std::move(*root), std::move(*reader), manifest->identity}};
#else
    static_cast<void>(parent);
    static_cast<void>(filename);
    static_cast<void>(source);
    static_cast<void>(manifest_limits);
    static_cast<void>(limits);
    static_cast<void>(plan);
    mark(diagnostic, ReferenceInputStage::Root);
    return fail(core::FailureCode::Unsupported);
#endif
}

core::Result<ReferenceInput>
ReferenceInput::from_memory(std::span<const std::byte> bytes, SourceId source,
                            std::optional<ReadRoot> root, const ReferenceInputLimits& limits,
                            core::ResourcePlan& plan, const std::stop_token& stop,
                            ReferenceInputDiagnostic* diagnostic) noexcept {
    clear(diagnostic);
    if (stop.stop_requested()) {
        return fail(core::FailureCode::Cancelled);
    }
    const auto document = read_document(bytes, source, limits.document, plan, stop, diagnostic);
    if (!document) {
        return std::unexpected(document.error());
    }
    auto schema = read_case(*document, limits.schema, plan, stop, diagnostic);
    if (!schema) {
        return std::unexpected(schema.error());
    }

    mark(diagnostic, ReferenceInputStage::Root);
#if defined(_WIN32) || defined(__linux__)
    if (root && !root->is_open()) {
        return invalid();
    }
#else
    if (root) {
        return fail(core::FailureCode::Unsupported);
    }
#endif
    if (!root && schema->data_.values()[0].records_FileReference.size() != 0) {
        return invalid();
    }

    if (cancelled(stop, diagnostic)) {
        return fail(core::FailureCode::Cancelled);
    }
    clear(diagnostic);
    return ReferenceInput{std::move(*schema), MemoryOrigin{std::move(root)}};
}

} // namespace uvlas::io::detail
