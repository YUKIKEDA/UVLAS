#pragma once

#include "reference_file_snapshot.hpp"
#include "reference_path.hpp"

#ifdef _WIN32
#include "windows_reference_root.hpp"
#elif defined(__linux__)
#include "linux_reference_root.hpp"
#endif

#include <filesystem>
#include <uvlas/io/case_schema.hpp>
#include <variant>

namespace uvlas::io::detail {

#ifdef _WIN32
using ReadRoot = WindowsReferenceRoot;
using ManifestReader = WindowsReferenceFile;
using ManifestOpenDiagnostic = WindowsReferenceOpenDiagnostic;
#elif defined(__linux__)
using ReadRoot = LinuxReferenceRoot;
using ManifestReader = LinuxReferenceFile;
using ManifestOpenDiagnostic = LinuxReferenceOpenReason;
#else
struct ReadRoot {};
struct ManifestReader {};
struct ManifestOpenDiagnostic {};
#endif

struct ReferenceInputLimits {
    InputLimits document;
    SchemaLimits schema;
};

struct ManifestLimits {
    std::size_t maximum_name_bytes;
    ReferenceFileLimits file;
};

enum class ReferenceInputStage : std::uint8_t {
    ManifestName,
    Root,
    ManifestOpen,
    ManifestRead,
    Document,
    Schema,
    Identity,
};

struct ReferenceInputDiagnostic {
    std::optional<ReferenceInputStage> stage;
    std::optional<ReferencePathDiagnostic> path;
    std::optional<ManifestOpenDiagnostic> open;
    std::optional<ReferenceFileDiagnostic> file;
    std::optional<InputDiagnostic> document;
    SchemaDiagnostic schema;

    void reset() noexcept;
};

// Borrows the schema and root from the input until it is moved or destroyed.
struct ReferenceInputView {
    const SchemaCase& schema;
    const ReadRoot* root = nullptr;
    std::optional<ReferenceFileId> manifest_identity;
};

// Owns the validated case together with its acquisition origin. Only the acquisition
// factories can pair a schema with a root or establish a native manifest identity.
class ReferenceInput {
public:
    ReferenceInput(ReferenceInput&&) noexcept = default;
    ReferenceInput& operator=(ReferenceInput&&) = delete;

    // The full location is specified as its native parent directory and one UTF-8
    // filename. Both the directory and the file are acquired by this call.
    [[nodiscard]] static core::Result<ReferenceInput>
    from_file(const std::filesystem::path& parent, std::string_view filename, SourceId source,
              ManifestLimits manifest_limits, const ReferenceInputLimits& limits,
              core::ResourcePlan& plan, const std::stop_token& stop,
              ReferenceInputDiagnostic* diagnostic) noexcept;

    // Consumes the caller's held directory, closing it if acquisition fails.
    // The root may be absent only when the checked case has no external references.
    [[nodiscard]] static core::Result<ReferenceInput>
    from_memory(std::span<const std::byte> bytes, SourceId source, std::optional<ReadRoot> root,
                const ReferenceInputLimits& limits, core::ResourcePlan& plan,
                const std::stop_token& stop, ReferenceInputDiagnostic* diagnostic) noexcept;

    [[nodiscard]] core::Result<ReferenceInputView> view() const& noexcept;
    [[nodiscard]] core::Result<ReferenceInputView> view() const&& = delete;

private:
    struct FileOrigin {
        ReadRoot root;
        ManifestReader reader;
        ReferenceFileId identity;
    };

    struct MemoryOrigin {
        std::optional<ReadRoot> root;
    };

    using Origin = std::variant<FileOrigin, MemoryOrigin>;

    ReferenceInput(SchemaCase schema, Origin origin) noexcept;

    SchemaCase schema_;
    Origin origin_;
};

} // namespace uvlas::io::detail
