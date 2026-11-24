#include "inspect_result.hpp"

#include "reference_input.hpp"
#include "result_documents.hpp"

#include <iostream>
#include <string_view>

#if defined(_WIN32)
#include <cstdio>
#include <fcntl.h>
#include <io.h>
#endif

namespace uvlas::app {
namespace {

using io::detail::FailedResultExit;

int reject(core::FailureCode reason, std::string_view name) {
    std::cerr << "inspect: ";
    if (reason == core::FailureCode::Cancelled) {
        std::cerr << "inspection interrupted\n";
        return static_cast<int>(FailedResultExit::Interrupted);
    }
    if (reason == core::FailureCode::ResourceExceeded) {
        std::cerr << "resource limit exceeded at " << name << '\n';
        return static_cast<int>(FailedResultExit::Resources);
    }
    std::cerr << "could not verify " << name << '\n';
    return static_cast<int>(FailedResultExit::Persistence);
}

#if defined(_WIN32) || defined(__linux__)

using namespace io::detail;

constexpr std::uint64_t kib = 1024;
constexpr std::uint64_t mib = 1024 * kib;
constexpr std::uint64_t metadata_cap = mib;
constexpr std::uint64_t resolved_cap = 64 * mib;
constexpr std::uint64_t payload_cap = 128 * mib;

[[nodiscard]] core::Result<ReferenceFileSnapshot>
read(const ReadRoot& root, std::string_view filename, ReferenceFileLimits limits,
     core::ResourcePlan& plan, const std::stop_token& stop,
     const io::SourceDigest* digest = nullptr) noexcept {
    const auto path = normalize_reference_path(filename, {filename.size(), filename.size(), 1},
                                               plan, stop, nullptr);
    if (!path) {
        return std::unexpected(path.error());
    }
    auto reader = root.open_file(*path, plan, stop, nullptr);
    if (!reader) {
        return std::unexpected(reader.error());
    }
    if (digest) {
        return read_reference_file(*reader, *digest, limits, plan, stop, nullptr);
    }
    return read_manifest_file(*reader, limits, plan, stop, nullptr);
}

#endif

} // namespace

int inspect_result(const std::filesystem::path& directory) {
#if defined(_WIN32) || defined(__linux__)
    const std::stop_token stop;
    const auto root = ReadRoot::open(directory, stop);
    if (!root) {
        return reject(root.error().code, "result directory");
    }
    core::ResourcePlan plan{payload_cap};
    constexpr ReferenceFileLimits metadata_file{metadata_cap, metadata_cap + 1, 4096, metadata_cap};
    const auto status_file = read(*root, "status.json", metadata_file, plan, stop);
    if (!status_file) {
        return reject(status_file.error().code, "status.json");
    }
    const auto provenance_file = read(*root, "provenance.json", metadata_file, plan, stop);
    if (!provenance_file) {
        return reject(provenance_file.error().code, "provenance.json");
    }
    const auto status = status_file->view();
    const auto provenance = provenance_file->view();
    if (!status || !provenance) {
        return reject(core::FailureCode::InvalidInput, "result metadata");
    }
    const io::InputLimits syntax{
        metadata_cap, 32,  8192,         8192,
        64 * kib,     128, metadata_cap, {io::ParserAllocationMode::EstimateAdmission, 32 * mib}};
    const auto checked =
        read_result_documents(status->bytes, provenance->bytes, syntax, plan, stop);
    if (!checked) {
        return reject(checked.error().code, "status.json and provenance.json");
    }
    if (const auto digest = checked->provenance().resolved_case) {
        const auto resolved =
            read(*root, "resolved-case.json", {resolved_cap, resolved_cap + 1, 8192, resolved_cap},
                 plan, stop, &*digest);
        if (!resolved) {
            return reject(resolved.error().code, "resolved-case.json");
        }
    }
    const std::string_view text{reinterpret_cast<const char*>(status->bytes.data()),
                                status->bytes.size()};
#if defined(_WIN32)
    if (_setmode(_fileno(stdout), _O_BINARY) == -1) {
        return reject(core::FailureCode::InvalidInput, "status output");
    }
#endif
    std::cout << text << std::flush;
    return std::cout ? 0 : reject(core::FailureCode::InvalidInput, "status output");
#else
    static_cast<void>(directory);
    return reject(core::FailureCode::Unsupported, "this platform's result directory reader");
#endif
}

} // namespace uvlas::app
