#include "reference_root_plan.hpp"

#include <algorithm>
#include <limits>
#include <utility>

namespace uvlas::io::detail {
namespace {

constexpr auto operation = core::Operation::PrepareModel;
using Reason = RootReferenceFailure;

[[nodiscard]] auto invalid() noexcept {
    return std::unexpected(core::Failure{core::FailureCode::InvalidInput, operation});
}

} // namespace

struct RootReferenceContext {
    const RootReferenceLimits& limits;
    core::ResourcePlan& plan;
    const std::stop_token& stop;
    RootReferenceDiagnostic* diagnostic = nullptr;
    std::optional<SchemaIdentity> parent{};
    std::optional<DeclaredReferenceView> current{};
    std::uint64_t work = 0;

    [[nodiscard]] std::unexpected<core::Failure>
    reject(core::Failure error, Reason reason, std::optional<std::uint64_t> required = {},
           std::optional<std::uint64_t> limit = {},
           std::optional<ReferencePathDiagnostic> path = {}) const noexcept {
        if (stop.stop_requested() || error.code == core::FailureCode::Cancelled) {
            return std::unexpected(core::Failure{core::FailureCode::Cancelled, operation});
        }
        if (diagnostic) {
            diagnostic->detail_ = RootReferenceDiagnosticDetail{
                reason,
                parent,
                current ? std::optional{current->source_id} : std::nullopt,
                current ? std::optional{current->declaration_node} : std::nullopt,
                path,
                required,
                limit};
            if (current && !current->declaration_pointer.empty() &&
                current->declaration_pointer.size() <= limits.diagnostic_bytes &&
                std::in_range<core::Index>(current->declaration_pointer.size())) {
                auto copy = plan.allocate<char>(
                    static_cast<core::Index>(current->declaration_pointer.size()), operation);
                if (copy) {
                    std::ranges::copy(current->declaration_pointer, copy->values().begin());
                    diagnostic->pointer_ = std::move(*copy);
                }
            }
            if (stop.stop_requested()) {
                diagnostic->reset();
                return std::unexpected(core::Failure{core::FailureCode::Cancelled, operation});
            }
        }
        return std::unexpected(error);
    }

    [[nodiscard]] std::unexpected<core::Failure>
    quota(Reason reason, std::optional<std::uint64_t> required,
          std::optional<std::uint64_t> limit = {}) const noexcept {
        return reject({core::FailureCode::ResourceExceeded, operation}, reason, required, limit);
    }

    [[nodiscard]] core::Result<void> charge(std::uint64_t count) noexcept {
        if (stop.stop_requested()) {
            return std::unexpected(core::Failure{core::FailureCode::Cancelled, operation});
        }
        if (count > limits.work - work) {
            const auto max = std::numeric_limits<std::uint64_t>::max();
            return quota(Reason::Work,
                         count <= max - work ? std::optional{work + count} : std::nullopt,
                         limits.work);
        }
        work += count;
        return {};
    }

    template <class T>
    [[nodiscard]] core::Result<core::PlannedBuffer<T>> allocate(core::Index count) const noexcept {
        if (stop.stop_requested()) {
            return std::unexpected(core::Failure{core::FailureCode::Cancelled, operation});
        }
        auto result = plan.allocate<T>(count, operation);
        if (stop.stop_requested()) {
            return std::unexpected(core::Failure{core::FailureCode::Cancelled, operation});
        }
        if (!result) {
            return reject(result.error(), Reason::Storage);
        }
        return result;
    }
};

const std::optional<RootReferenceDiagnosticDetail>&
RootReferenceDiagnostic::detail() const& noexcept {
    return detail_;
}

std::optional<std::string_view> RootReferenceDiagnostic::pointer() const& noexcept {
    return pointer_.size()
               ? std::optional{std::string_view{pointer_.values().data(), pointer_.size()}}
               : std::nullopt;
}

void RootReferenceDiagnostic::reset() noexcept {
    detail_.reset();
    pointer_ = {};
}

PreparedRootReferences::PreparedRootReferences(DeclaredReferences declarations,
                                               core::PlannedBuffer<core::Index> edges,
                                               core::PlannedBuffer<std::optional<Source>> sources,
                                               RootReferenceReport report) noexcept
    : declarations_(std::move(declarations)),
      edges_(std::move(edges)),
      sources_(std::move(sources)),
      report_(report) {}

core::Result<SchemaIdentity> PreparedRootReferences::identity() const noexcept {
    return declarations_.identity();
}

core::Result<RootReferenceReport> PreparedRootReferences::report() const noexcept {
    if (!identity()) {
        return invalid();
    }
    return report_;
}

core::Result<RootReferenceDeclarationView>
PreparedRootReferences::declaration(core::Index index) const& noexcept {
    const auto original = declarations_.entry(index);
    if (!original || index < 0 || std::cmp_greater_equal(index, edges_.size())) {
        return invalid();
    }
    return RootReferenceDeclarationView{*original,
                                        edges_.values()[static_cast<std::size_t>(index)]};
}

core::Result<RootReferenceSourceView>
PreparedRootReferences::source(core::Index index) const& noexcept {
    if (index < 0 || std::cmp_greater_equal(index, sources_.size())) {
        return invalid();
    }
    const auto& slot = sources_.values()[static_cast<std::size_t>(index)];
    if (!slot) {
        return invalid();
    }
    return RootReferenceSourceView{slot->first_declaration, slot->path};
}

core::Result<PreparedRootReferences>
prepare_root_references(DeclaredReferences&& declarations, RootReferenceLimits limits,
                        core::ResourcePlan& plan, const std::stop_token& stop,
                        RootReferenceDiagnostic* diagnostic) noexcept {
    if (diagnostic) {
        diagnostic->reset();
    }
    if (stop.stop_requested()) {
        return std::unexpected(core::Failure{core::FailureCode::Cancelled, operation});
    }
    const auto allocation_start = plan.allocated_bytes();
    RootReferenceContext context{limits, plan, stop, diagnostic};
    const auto identity = declarations.identity();
    const auto count = declarations.size();
    if (!identity || !count) {
        return context.reject(invalid().error(), Reason::InvalidOwner);
    }
    context.parent = *identity;

    if (*count > limits.declarations) {
        return context.quota(Reason::Declarations, *count, limits.declarations);
    }
    if (*count > limits.root_edges) {
        return context.quota(Reason::RootEdges, *count, limits.root_edges);
    }
    const auto capacity = std::min(*count, limits.unique_sources);
    constexpr auto max_index = static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max());
    if (*count > max_index || capacity > (max_index - 1) / 2) {
        return context.quota(Reason::Storage, std::nullopt);
    }
    if (*count != 0 && capacity == 0) {
        return context.quota(Reason::UniqueSources, 1, limits.unique_sources);
    }
    const auto buckets = capacity == 0 ? 0 : 2 * capacity + 1;
    for (const auto units : {*count, capacity, buckets}) {
        if (const auto charged = context.charge(units); !charged) {
            return std::unexpected(charged.error());
        }
    }

    auto edges = context.allocate<core::Index>(static_cast<core::Index>(*count));
    if (!edges) {
        return std::unexpected(edges.error());
    }
    auto paths = context.allocate<std::optional<PreparedRootReferences::Source>>(
        static_cast<core::Index>(capacity));
    if (!paths) {
        return std::unexpected(paths.error());
    }
    // An empty bucket is zero. Entries hold one plus the first-seen source ordinal.
    auto table = context.allocate<std::uint64_t>(static_cast<core::Index>(buckets));
    if (!table) {
        return std::unexpected(table.error());
    }
    std::uint64_t unique = 0;
    std::uint64_t normalized_bytes = 0;
    for (std::uint64_t i = 0; i < *count; ++i) {
        if (const auto charged = context.charge(1); !charged) {
            return std::unexpected(charged.error());
        }
        const auto occurrence = declarations.entry(static_cast<core::Index>(i));
        if (!occurrence) {
            return context.reject(occurrence.error(), Reason::InvalidOwner);
        }
        context.current = *occurrence;
        std::uint64_t slot = occurrence->source_id.value() % buckets;
        for (;;) {
            if (const auto charged = context.charge(1); !charged) {
                return std::unexpected(charged.error());
            }
            auto& entry = table->values()[static_cast<std::size_t>(slot)];
            if (entry == 0) {
                if (unique == capacity) {
                    return context.quota(Reason::UniqueSources, unique + 1, limits.unique_sources);
                }
                constexpr std::uint64_t passes = 16;
                const auto length = occurrence->original_path.size();
                if (length > std::numeric_limits<std::uint64_t>::max() / passes) {
                    return context.quota(Reason::Work, std::nullopt, limits.work);
                }
                // The existing normalizer does bounded forward and backward scans and one copy.
                // Charge a conservative byte-pass envelope before calling it.
                if (const auto charged = context.charge(passes * length); !charged) {
                    return std::unexpected(charged.error());
                }
                const auto remaining = limits.normalized_path_bytes - normalized_bytes;
                std::optional<ReferencePathDiagnostic> path_diagnostic;
                auto normalized = normalize_reference_path(
                    occurrence->original_path,
                    {length,
                     static_cast<std::size_t>(std::min<std::uint64_t>(
                         remaining, std::numeric_limits<std::size_t>::max())),
                     limits.path_components},
                    plan, stop, &path_diagnostic);
                if (!normalized) {
                    return context.reject(normalized.error(), Reason::Path, {}, {},
                                          path_diagnostic);
                }
                normalized_bytes += normalized->text().size();
                paths->values()[static_cast<std::size_t>(unique)].emplace(
                    PreparedRootReferences::Source{static_cast<core::Index>(i),
                                                   std::move(*normalized)});
                entry = ++unique;
                edges->values()[static_cast<std::size_t>(i)] = static_cast<core::Index>(entry - 1);
                break;
            }
            const auto& prior = paths->values()[static_cast<std::size_t>(entry - 1)];
            if (!prior) {
                return context.reject(invalid().error(), Reason::InvalidOwner);
            }
            const auto first = declarations.entry(prior->first_declaration);
            if (!first) {
                return context.reject(first.error(), Reason::InvalidOwner);
            }
            if (first->source_id == occurrence->source_id) {
                edges->values()[static_cast<std::size_t>(i)] = static_cast<core::Index>(entry - 1);
                break;
            }
            // At most half of the buckets can be occupied, so a probe always terminates.
            slot = slot + 1 == buckets ? 0 : slot + 1;
        }
        context.current.reset();
    }
    if (stop.stop_requested()) {
        return std::unexpected(core::Failure{core::FailureCode::Cancelled, operation});
    }
    return PreparedRootReferences{std::move(declarations),
                                  std::move(*edges),
                                  std::move(*paths),
                                  {*count, unique, normalized_bytes, context.work,
                                   plan.allocated_bytes() - allocation_start}};
}

} // namespace uvlas::io::detail
