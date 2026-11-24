#include "validate_case.hpp"

#include "reference_resolution.hpp"
#include "result_documents.hpp"

#include <cstdint>
#include <iostream>
#include <stop_token>
#include <string>
#include <string_view>
#include <utility>

namespace uvlas::app {
namespace {

using namespace io::detail;

// Host admission limits apply before the manifest's own quotas can be trusted.
constexpr std::uint64_t kib = 1024;
constexpr std::uint64_t mib = 1024 * kib;
constexpr std::uint64_t payload_cap = 256 * mib;
constexpr std::uint64_t manifest_cap = 8 * mib;
constexpr std::uint64_t reference_cap = 32 * mib;
constexpr std::uint64_t reference_total_cap = 64 * mib;
constexpr std::uint64_t maximum_references = 1024;
constexpr std::uint64_t maximum_work = 64 * mib;
constexpr std::uint64_t maximum_nodes = 256 * kib;
constexpr std::size_t maximum_path = 4096;

ReferenceInputLimits input_limits() noexcept {
    return {
        {manifest_cap,
         128,
         maximum_nodes,
         maximum_nodes,
         64 * kib,
         1024,
         manifest_cap,
         {io::ParserAllocationMode::EstimateAdmission, payload_cap}},
        io::SchemaLimits{{
            manifest_cap,  128,           maximum_nodes, maximum_nodes, 64 * kib,
            1024,          manifest_cap,  64 * mib,      maximum_nodes, maximum_nodes,
            maximum_nodes, maximum_nodes, maximum_nodes, maximum_nodes, maximum_nodes,
            128,           1024,          maximum_nodes, maximum_path,  64 * mib,
            maximum_work,  4096,          4096,
        }},
    };
}

int report(FailedResultExit exit, std::string_view message, std::string_view location = {}) {
    std::cerr << "validate: " << message;
    if (!location.empty()) {
        std::cerr << " at " << location;
    }
    std::cerr << '\n';
    return static_cast<int>(exit);
}

int reject(const core::Failure& failure, std::string_view stage, std::string_view location = {}) {
    if (failure.code == core::FailureCode::ResourceExceeded) {
        return report(FailedResultExit::Resources, "resource limit exceeded", stage);
    }
    if (failure.code == core::FailureCode::Cancelled) {
        return report(FailedResultExit::Interrupted, "validation interrupted");
    }
    return report(FailedResultExit::Input, stage, location);
}

int reject_manifest(const core::Failure& failure, const ReferenceInputDiagnostic& diagnostic) {
    if (diagnostic.stage == ReferenceInputStage::Document) {
        if (const auto byte = diagnostic.document ? diagnostic.document->byte_offset : std::nullopt;
            byte) {
            return reject(failure, "JSON document at byte " + std::to_string(*byte));
        }
        return reject(failure, "JSON document");
    }
    if (diagnostic.stage == ReferenceInputStage::Schema) {
        return reject(failure, "case schema", diagnostic.schema.path().value_or(""));
    }
    return reject(failure, "manifest acquisition");
}

int validate_references(const ReferenceInput& input, const ReferenceInputView& view,
                        core::ResourcePlan& plan, const std::stop_token& stop) {
    const DeclarationLimits declaration_limits{{maximum_nodes, 128, maximum_references,
                                                manifest_cap, manifest_cap, 16 * mib, 16 * mib,
                                                32 * mib, maximum_work, maximum_path}};
    DeclarationDiagnostic declaration_diagnostic;
    auto declarations = collect_case_file_references(view.schema, declaration_limits, plan, stop,
                                                     &declaration_diagnostic);
    if (!declarations) {
        return reject(declarations.error(), "reference declarations",
                      declaration_diagnostic.path().value_or(""));
    }
    RootReferenceDiagnostic plan_diagnostic;
    auto prepared =
        prepare_root_references(std::move(*declarations),
                                {maximum_references, maximum_references, maximum_references,
                                 manifest_cap, 128, maximum_work, maximum_path},
                                plan, stop, &plan_diagnostic);
    if (!prepared) {
        return reject(prepared.error(), "reference paths", plan_diagnostic.pointer().value_or(""));
    }
    RootSnapshotDiagnostic snapshot_diagnostic;
    LeafReferenceDiagnostic diagnostic;
    const auto closed = resolve_leaf_references(
        input, std::move(*prepared),
        {{{reference_cap, reference_cap + 1, 8192, reference_cap},
          maximum_references,
          maximum_path},
         maximum_references,
         reference_total_cap,
         reference_total_cap + maximum_references,
         32 * kib,
         reference_total_cap,
         maximum_references * maximum_references},
        {maximum_references, maximum_references, 1, maximum_work, maximum_path},
        {reference_cap, maximum_nodes, maximum_nodes, 4 * maximum_nodes, maximum_nodes,
         maximum_nodes, maximum_work},
        plan, stop, &snapshot_diagnostic, &diagnostic);
    if (!closed) {
        if (diagnostic.reason == LeafReferenceFailure::UnsupportedRole) {
            return reject(closed.error(), "external reference dependency validation is unsupported",
                          diagnostic.pointer());
        }
        if (diagnostic.reason == LeafReferenceFailure::Format) {
            const auto message =
                diagnostic.format_offset
                    ? "Gmsh mesh format at byte " + std::to_string(*diagnostic.format_offset)
                    : "Gmsh mesh format";
            return reject(closed.error(), message, diagnostic.pointer());
        }
        if (diagnostic.reason == LeafReferenceFailure::Graph) {
            return reject(closed.error(), "external reference dependency graph",
                          diagnostic.pointer());
        }
        const auto message = snapshot_diagnostic.read && snapshot_diagnostic.read->reason ==
                                                             ReferenceFileFailure::DigestMismatch
                                 ? "reference SHA-256 mismatch"
                                 : "reference acquisition";
        return reject(closed.error(), message, snapshot_diagnostic.pointer());
    }
    std::cout << "Case validation succeeded.\n";
    return 0;
}

} // namespace

int validate_case(const std::filesystem::path& path) {
    auto parent = path.parent_path();
    if (parent.empty()) {
        parent = ".";
    }
    const auto filename = path.filename().u8string();
    const std::string_view name{reinterpret_cast<const char*>(filename.data()), filename.size()};
    core::ResourcePlan plan{payload_cap};
    const std::stop_token stop;
    ReferenceInputDiagnostic diagnostic;
    const auto input = ReferenceInput::from_file(
        parent, name, io::SourceId{0},
        {maximum_path, {manifest_cap, manifest_cap + 1, 4096, manifest_cap}}, input_limits(), plan,
        stop, &diagnostic);
    if (!input) {
        return reject_manifest(input.error(), diagnostic);
    }
    const auto view = input->view();
    if (!view) {
        return reject(view.error(), "manifest identity");
    }
    return validate_references(*input, *view, plan, stop);
}

} // namespace uvlas::app
