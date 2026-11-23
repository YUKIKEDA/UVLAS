#pragma once

#include <optional>
#include <stop_token>
#include <string_view>
#include <uvlas/core/state.hpp>
#include <uvlas/io/input_document.hpp>

namespace uvlas::io::detail {

// The caller chooses the termination category in its operation context. The
// original Failure keeps its independent cause and operation.
enum class FailedResultExit : std::uint8_t {
    Input = 2,
    Geometry = 3,
    Numerical = 4,
    Resources = 5,
    Persistence = 6,
    Interrupted = 130
};

struct AcceptedResultIdentity {
    core::AcceptedStateId id;
    core::ModelGeneration model_generation;
    core::NumericalGeneration numerical_generation;
    core::ModelFingerprint model;
};

class ResultStatus {
public:
    [[nodiscard]] static core::Result<ResultStatus>
    completed(const core::AcceptedState* last = nullptr) noexcept;
    [[nodiscard]] static core::Result<ResultStatus>
    failed(core::Failure, FailedResultExit, const core::AcceptedState* last = nullptr) noexcept;

    [[nodiscard]] std::uint8_t exit_code() const noexcept {
        return exit_code_;
    }
    [[nodiscard]] std::optional<core::Failure> failure() const noexcept {
        return failure_;
    }
    [[nodiscard]] std::optional<AcceptedResultIdentity> accepted() const noexcept {
        return accepted_;
    }

private:
    ResultStatus(std::uint8_t exit_code, std::optional<core::Failure> failure,
                 std::optional<AcceptedResultIdentity> accepted) noexcept;

    std::uint8_t exit_code_;
    std::optional<core::Failure> failure_;
    std::optional<AcceptedResultIdentity> accepted_;
};

// Identities supplied by the caller's source and resolved-case producers.
// No file is opened, read, or implicitly declared resolved by this projection.
struct ResultProvenance {
    std::optional<SourceIdentity> original_case;
    std::optional<SourceDigest> resolved_case;
};

class ResultDocuments {
public:
    ResultDocuments(ResultDocuments&&) noexcept = default;
    ResultDocuments& operator=(ResultDocuments&&) noexcept = default;

    [[nodiscard]] std::string_view status() const& noexcept;
    [[nodiscard]] std::string_view status() const&& = delete;
    [[nodiscard]] std::string_view provenance() const& noexcept;
    [[nodiscard]] std::string_view provenance() const&& = delete;

private:
    friend core::Result<ResultDocuments> write_result_documents(const ResultStatus&,
                                                                const ResultProvenance&,
                                                                std::uint64_t, core::ResourcePlan&,
                                                                const std::stop_token&) noexcept;

    ResultDocuments(core::PlannedBuffer<char> bytes, std::size_t status_bytes) noexcept;

    core::PlannedBuffer<char> bytes_;
    std::size_t status_bytes_;
};

// Private UTF-8 projection for status.json and provenance.json, with SHA-256
// linking the exact provenance bytes. Charges the combined owned output once.
// Serialization can fail independently of the original operation outcome.
[[nodiscard]] core::Result<ResultDocuments>
write_result_documents(const ResultStatus&, const ResultProvenance&,
                       std::uint64_t maximum_output_bytes, core::ResourcePlan&,
                       const std::stop_token&) noexcept;

} // namespace uvlas::io::detail
