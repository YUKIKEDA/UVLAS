#pragma once

#include "result_documents.hpp"

#include <uvlas/io/si_case.hpp>

namespace uvlas::io::detail {

struct ResolvedCaseDocumentLimits {
    std::uint64_t maximum_quantities;
    std::uint64_t maximum_output_bytes;
    // Work covers both output passes and hashing the completed bytes.
    std::uint64_t maximum_work;
};

class ResolvedCaseDocument {
public:
    ResolvedCaseDocument(ResolvedCaseDocument&&) noexcept = default;
    ResolvedCaseDocument& operator=(ResolvedCaseDocument&&) noexcept = default;

    [[nodiscard]] std::string_view text() const& noexcept;
    [[nodiscard]] std::string_view text() const&& = delete;
    [[nodiscard]] std::optional<ResultProvenance> provenance() const noexcept;

private:
    friend core::Result<ResolvedCaseDocument>
    write_resolved_case_document(const SiCase&, ResolvedCaseDocumentLimits, core::ResourcePlan&,
                                 const std::stop_token&) noexcept;

    ResolvedCaseDocument(core::PlannedBuffer<char> text, SourceIdentity source,
                         SourceDigest digest) noexcept;

    core::PlannedBuffer<char> text_;
    SourceIdentity source_;
    SourceDigest digest_;
};

struct ResolvedResultDocuments {
    ResolvedCaseDocument resolved_case;
    ResultDocuments result;
};

// Private report of the complete checked input, SI ledger and actual structural defaults.
// Numeric output contains the existing enclosures, including both retained words. Non-SI
// preparation and a public persisted schema are separate from this byte projection.
[[nodiscard]] core::Result<ResolvedCaseDocument>
write_resolved_case_document(const SiCase&, ResolvedCaseDocumentLimits, core::ResourcePlan&,
                             const std::stop_token&) noexcept;

// The provenance digest is derived from the exact resolved output by this call. No
// partial result is returned if any writer fails; successful Plan allocations remain charged.
[[nodiscard]] core::Result<ResolvedResultDocuments>
write_resolved_result_documents(const SiCase&, const ResultStatus&, ResolvedCaseDocumentLimits,
                                std::uint64_t maximum_result_bytes, core::ResourcePlan&,
                                const std::stop_token&) noexcept;

} // namespace uvlas::io::detail
