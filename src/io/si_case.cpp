#include "si_case_detail.hpp"

#include <charconv>
#include <new>

namespace uvlas::io {
namespace {
core::Failure invalid_owner() noexcept {
    return {core::FailureCode::InvalidInput, core::Operation::PrepareModel};
}
} // namespace
core::Result<std::uint64_t> SiLimits::get(SiQuota quota) const noexcept {
    const auto i = static_cast<std::size_t>(quota);
    if (i >= quotas_.size())
        return std::unexpected(invalid_owner());
    return quotas_[i];
}
bool SiQuantityKey::operator==(const SiQuantityKey& other) const noexcept {
    return owner_ == other.owner_ && detail::same_identity(source_, other.source_) &&
           node_ == other.node_ && type_ == other.type_ && role_ == other.role_;
}
SiCatalog::SiCatalog(core::PlannedBuffer<detail::SiCatalogData> p) noexcept : data_(std::move(p)) {}
SiCatalog::SiCatalog(SiCatalog&&) noexcept = default;
SiCatalog& SiCatalog::operator=(SiCatalog&&) noexcept = default;
SiCatalog::~SiCatalog() = default;
core::Result<std::uint64_t> SiCatalog::size() const noexcept {
    if (!data_.size())
        return std::unexpected(invalid_owner());
    return data_.values()[0].entries.size();
}
core::Result<SiCatalogEntry> SiCatalog::entry(std::uint64_t row) const& noexcept {
    if (!data_.size() || row >= data_.values()[0].entries.size())
        return std::unexpected(invalid_owner());
    return detail::si_catalog_entry(data_.values()[0], row);
}
core::Result<SiMemoryReport> SiCatalog::memory_report() const noexcept {
    if (!data_.size())
        return std::unexpected(invalid_owner());
    return data_.values()[0].report;
}
SiRequest::SiRequest(core::PlannedBuffer<detail::SiRequestData> p) noexcept : data_(std::move(p)) {}
SiRequest::SiRequest(SiRequest&&) noexcept = default;
SiRequest& SiRequest::operator=(SiRequest&&) noexcept = default;
SiRequest::~SiRequest() = default;
core::Result<SiMemoryReport> SiRequest::memory_report() const noexcept {
    if (!data_.size())
        return std::unexpected(invalid_owner());
    return data_.values()[0].report;
}
SiCase::SiCase(core::PlannedBuffer<detail::SiData> p) noexcept : data_(std::move(p)) {}
SiCase::SiCase(SiCase&&) noexcept = default;
SiCase& SiCase::operator=(SiCase&&) noexcept = default;
SiCase::~SiCase() = default;
core::Result<std::reference_wrapper<const SchemaCase>> SiCase::source_case() const& noexcept {
    if (!data_.size() || !data_.values()[0].source)
        return std::unexpected(invalid_owner());
    return std::cref(*data_.values()[0].source);
}
core::Result<SiResolvedDefaults> SiCase::defaults() const noexcept {
    if (!data_.size())
        return std::unexpected(invalid_owner());
    return data_.values()[0].defaults;
}
core::Result<SiMemoryReport> SiCase::memory_report() const noexcept {
    if (!data_.size())
        return std::unexpected(invalid_owner());
    return data_.values()[0].report;
}
core::Result<SiQuantityView> SiCase::quantity_impl(const detail::SchemaData* owner,
                                                   std::uint64_t node) const noexcept {
    if (!data_.size() || !data_.values()[0].source)
        return std::unexpected(invalid_owner());
    const auto& data = data_.values()[0];
    if (!data.source.has_value())
        return std::unexpected(invalid_owner());
    if (owner != detail::SiFactory::storage(*data.source) || node >= data.node_to_row.size())
        return std::unexpected(invalid_owner());
    const auto row = data.node_to_row.values()[static_cast<std::size_t>(node)];
    if (row >= data.entries.size())
        return std::unexpected(invalid_owner());
    const auto& entry = data.entries.values()[static_cast<std::size_t>(row)];
    const auto identity = data.source->identity();
    if (!identity)
        return std::unexpected(identity.error());
    return SiQuantityView{
        detail::SiFactory::key(
            data.key_owner, *identity, node,
            static_cast<std::uint16_t>(owner->nodes.values()[static_cast<std::size_t>(node)].type),
            entry.record),
        detail::si_recipe(owner->snapshot.values()[0], owner->nodes.values(),
                          owner->digits.values(), node, data.units),
        data.converted.values()[static_cast<std::size_t>(row)], row,
        detail::si_requirement_view(data.requirements.values()[static_cast<std::size_t>(row)],
                                    data.request_text.values())};
}
std::optional<std::string_view> SiDiagnostic::path() const& noexcept {
    if (!detail_ || detail_->path_state != SchemaPathState::Complete)
        return std::nullopt;
    return std::string_view{path_.values().data(), path_.size()};
}
void SiDiagnostic::reset() noexcept {
    detail_.reset();
    path_ = {};
    requirement_.reset();
    requirement_text_ = {};
}
SiDiagnostic::SiDiagnostic(SiDiagnostic&& other) noexcept {
    *this = std::move(other);
}
SiDiagnostic& SiDiagnostic::operator=(SiDiagnostic&& other) noexcept {
    if (this != &other) {
        detail_ = std::move(other.detail_);
        path_ = std::move(other.path_);
        requirement_ = other.requirement_;
        requirement_text_ = std::move(other.requirement_text_);
        other.reset();
    }
    return *this;
}
std::optional<SiRequirementView> SiDiagnostic::requirement() const& noexcept {
    return requirement_;
}

namespace detail {
SiContext::SiContext(const SiLimits& configured, core::ResourcePlan& plan, SiDiagnostic* diagnostic)
    : limits(configured),
      resources(plan),
      slot(diagnostic),
      work{0, *configured.get(SiQuota::Work), configured.stop()},
      start_charged(plan.allocated_bytes()) {
    if (slot)
        slot->reset();
}
std::uint64_t SiContext::limit(SiQuota quota) const {
    return *limits.get(quota);
}
std::uint64_t SiContext::add(std::uint64_t a, std::uint64_t b) {
    if (b > std::numeric_limits<std::uint64_t>::max() - a)
        reject(SiReason::Overflow, core::FailureCode::ResourceExceeded);
    return a + b;
}
std::uint64_t SiContext::multiply(std::uint64_t a, std::uint64_t b) {
    if (a && b > std::numeric_limits<std::uint64_t>::max() / a)
        reject(SiReason::Overflow, core::FailureCode::ResourceExceeded);
    return a * b;
}
void SiContext::require(SiQuota quota, std::uint64_t count) {
    if (count > limit(quota))
        reject(SiReason::QuotaExceeded, core::FailureCode::ResourceExceeded, quota, count);
}
void SiContext::admit(std::uint64_t p, std::uint64_t s, std::uint64_t w, std::uint64_t path,
                      std::optional<std::uint64_t> allocated_scratch,
                      std::uint64_t cumulative_extra, std::uint64_t diagnostic_request_bytes) {
    planned_owned = p;
    planned_scratch = s;
    planned_work = add(add(work.used, w), multiply(add(path, diagnostic_request_bytes), 2));
    admission = add(add(add(p, allocated_scratch.value_or(s)), cumulative_extra),
                    add(path, diagnostic_request_bytes));
    require(SiQuota::OwnedBytes, p);
    require(SiQuota::ScratchBytes, s);
    require(SiQuota::Work, planned_work);
    require(SiQuota::DiagnosticPathBytes, path);
    require(SiQuota::CumulativeBytes, admission);
    if (admission > resources.remaining_bytes())
        reject(SiReason::QuotaExceeded, core::FailureCode::ResourceExceeded,
               SiQuota::CumulativeBytes, admission);
}
namespace {
std::string_view key_text(const InputDocumentData& source, std::uint64_t member) {
    const auto& range = source.members.values()[static_cast<std::size_t>(member)].key;
    return {source.text.values().data() + range.offset, range.size};
}
std::uint64_t segment_bytes(SiContext& c, const SchemaNode& node) {
    std::uint64_t bytes = 1;
    if (node.parent_member) {
        for (const char ch : key_text(*c.snapshot, *node.parent_member)) {
            c.work.step();
            bytes = c.add(bytes, (ch == '~' || ch == '/') ? 2 : 1);
        }
    } else if (node.parent_array_index) {
        auto value = *node.parent_array_index;
        do {
            c.work.step();
            ++bytes;
            value /= 10;
        } while (value);
    } else
        throw invalid_owner();
    return bytes;
}
void capture_requirement(SiContext& c, SiDiagnostic& destination,
                         SiDiagnosticDetail& fixed) noexcept {
    if (!c.requirement)
        return;
    try {
        const auto& wanted = *c.requirement;
        const auto count = [&c](const SiExactInput& value) {
            if (const auto* decimal = std::get_if<SiDecimalInput>(&value))
                return static_cast<std::uint64_t>(decimal->token.size());
            const auto& rational = std::get<SiRationalInput>(value);
            return c.add(rational.numerator.size(), rational.denominator.size());
        };
        const auto bytes =
            c.add(c.add(count(wanted.scale.value), count(wanted.absolute)), count(wanted.relative));
        const auto spent = c.resources.allocated_bytes() - c.start_charged;
        if (bytes > c.limit(SiQuota::RequestTokenBytes) ||
            spent > c.limit(SiQuota::CumulativeBytes) ||
            bytes > c.limit(SiQuota::CumulativeBytes) - spent ||
            bytes > c.resources.remaining_bytes() ||
            bytes > static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max())) {
            fixed.request_state = SiRequestState::Quota;
            return;
        }
        c.work.step(c.multiply(bytes, 2));
        auto storage = c.resources.allocate<char>(static_cast<core::Index>(bytes),
                                                  core::Operation::PrepareModel);
        if (!storage) {
            fixed.request_state = SiRequestState::Allocation;
            return;
        }
        std::size_t offset = 0;
        const auto save = [&](std::string_view value) {
            auto* target = storage->values().data() + offset;
            std::ranges::copy(value, target);
            offset += value.size();
            return std::string_view{target, value.size()};
        };
        const auto copy = [&](const SiExactInput& value) -> SiExactInput {
            if (const auto* decimal = std::get_if<SiDecimalInput>(&value))
                return SiDecimalInput{save(decimal->token)};
            const auto& rational = std::get<SiRationalInput>(value);
            return SiRationalInput{save(rational.numerator), save(rational.denominator)};
        };
        const SiRequirementView view{{copy(wanted.scale.value), wanted.scale.dimensions},
                                     copy(wanted.absolute),
                                     copy(wanted.relative)};
        SiFactory::diagnostic_requirement(destination, std::move(*storage), view);
        fixed.request_state = SiRequestState::Complete;
    } catch (const std::bad_alloc&) {
        fixed.request_state = SiRequestState::Allocation;
    } catch (...) {
        fixed.request_state = SiRequestState::Quota;
    }
}
} // namespace
void SiContext::finish() {
    if (owned > planned_owned)
        reject(SiReason::QuotaExceeded, core::FailureCode::ResourceExceeded, SiQuota::OwnedBytes,
               owned);
    const auto peak = add(add(scratch, automatic_scratch), exact_scratch_peak);
    if (peak > planned_scratch)
        reject(SiReason::QuotaExceeded, core::FailureCode::ResourceExceeded, SiQuota::ScratchBytes,
               peak);
    if (work.used > planned_work)
        reject(SiReason::QuotaExceeded, core::FailureCode::ResourceExceeded, SiQuota::Work,
               work.used);
    if (resources.allocated_bytes() - start_charged > admission)
        reject(SiReason::QuotaExceeded, core::FailureCode::ResourceExceeded,
               SiQuota::CumulativeBytes, resources.allocated_bytes() - start_charged);
}
[[noreturn]] void SiContext::reject(SiReason reason, core::FailureCode code,
                                    std::optional<SiQuota> quota,
                                    std::optional<std::uint64_t> required) {
    if (slot) {
        auto* const destination = slot;
        SiDiagnosticDetail fixed{reason,
                                 source,
                                 key,
                                 quota,
                                 required,
                                 quota ? std::optional{limit(*quota)} : std::nullopt,
                                 dimensions,
                                 achieved,
                                 work.used,
                                 SiRequestState::NotAvailable,
                                 SchemaPathState::NotBuilt};
        SiFactory::diagnostic(*slot, fixed);
        // A diagnostic cannot recursively replace the failure it is trying to explain.
        slot = nullptr;
        bool complete = false;
        capture_requirement(*this, *destination, fixed);
        try {
            if (snapshot && current && *current < nodes.size()) {
                std::uint64_t count = 0, node = *current;
                while (nodes[static_cast<std::size_t>(node)].parent) {
                    work.step();
                    count = add(count, segment_bytes(*this, nodes[static_cast<std::size_t>(node)]));
                    node = *nodes[static_cast<std::size_t>(node)].parent;
                }
                if (count > limit(SiQuota::DiagnosticPathBytes))
                    fixed.path_state = SchemaPathState::PathQuota;
                else if (resources.allocated_bytes() - start_charged <=
                             limit(SiQuota::CumulativeBytes) &&
                         count <= limit(SiQuota::CumulativeBytes) -
                                      (resources.allocated_bytes() - start_charged) &&
                         count <=
                             static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max()) &&
                         count <= resources.remaining_bytes()) {
                    auto path = resources.allocate<char>(static_cast<core::Index>(count),
                                                         core::Operation::PrepareModel);
                    if (!path)
                        fixed.path_state = SchemaPathState::Allocation;
                    else {
                        std::uint64_t offset = count;
                        node = *current;
                        while (nodes[static_cast<std::size_t>(node)].parent) {
                            const auto& info = nodes[static_cast<std::size_t>(node)];
                            if (info.parent_member) {
                                const auto text = key_text(*snapshot, *info.parent_member);
                                for (std::size_t i = text.size(); i > 0; --i) {
                                    work.step();
                                    const char ch = text[i - 1];
                                    if (ch == '~' || ch == '/') {
                                        path->values()[--offset] = ch == '~' ? '0' : '1';
                                        path->values()[--offset] = '~';
                                    } else
                                        path->values()[--offset] = ch;
                                }
                            } else {
                                auto value = *info.parent_array_index;
                                do {
                                    work.step();
                                    path->values()[--offset] = static_cast<char>('0' + value % 10);
                                    value /= 10;
                                } while (value);
                            }
                            path->values()[--offset] = '/';
                            node = *info.parent;
                        }
                        fixed.path_state = SchemaPathState::Complete;
                        fixed.work = work.used;
                        SiFactory::diagnostic(*destination, fixed, std::move(*path));
                        complete = true;
                    }
                }
            }
        } catch (const std::bad_alloc&) {
            fixed.path_state = SchemaPathState::Allocation;
        } catch (...) {
            fixed.path_state = SchemaPathState::NotBuilt;
        }
        slot = destination;
        fixed.work = work.used;
        if (!complete)
            SiFactory::diagnostic(*slot, fixed);
    }
    throw core::Failure{code, core::Operation::PrepareModel};
}
} // namespace detail
} // namespace uvlas::io
