#pragma once

#include "case_schema_storage.hpp"
#include "si_interval.hpp"
#include "si_program.hpp"

#include <limits>
#include <uvlas/io/si_case.hpp>

namespace uvlas::io::detail {
struct SiKeyOwner {};
struct SiEntryStorage {
    std::uint64_t node = 0;
    std::uint64_t parent_record = 0;
    std::uint16_t record = 0;
    std::uint16_t field = 0;
};
struct SiCatalogData {
    std::shared_ptr<const SiKeyOwner> key_owner;
    std::optional<SchemaIdentity> identity;
    SiUnits units{SiLengthUnit::Metre, SiMassUnit::Kilogram, SiTimeUnit::Second,
                  SiAngleUnit::Radian};
    core::PlannedBuffer<InputDocumentData> snapshot;
    core::PlannedBuffer<SchemaNode> nodes;
    core::PlannedBuffer<char> digits;
    core::PlannedBuffer<SiEntryStorage> entries;
    std::uint64_t maximum_path_bytes = 0;
    SiMemoryReport report{};
};
struct SiTextRange {
    std::uint64_t offset = 0, count = 0;
};
struct SiExactStorage {
    struct Decimal {
        SiTextRange token;
    };
    struct Rational {
        SiTextRange numerator, denominator;
    };
    std::variant<Decimal, Rational> value;
};
struct SiRequirementStorage {
    std::uint64_t source_node = 0;
    SchemaDimension dimensions{{0, 0, 0, 0}};
    SiExactStorage scale, absolute, relative;
};
struct SiRequestData {
    std::shared_ptr<const SiKeyOwner> key_owner;
    std::optional<SchemaIdentity> identity;
    core::PlannedBuffer<SiRequirementStorage> rows;
    core::PlannedBuffer<char> text;
    SiMemoryReport report{};
};
struct SiData {
    std::shared_ptr<const SiKeyOwner> key_owner;
    std::optional<SchemaCase> source;
    SiUnits units{SiLengthUnit::Metre, SiMassUnit::Kilogram, SiTimeUnit::Second,
                  SiAngleUnit::Radian};
    core::PlannedBuffer<SiEntryStorage> entries;
    core::PlannedBuffer<SiRequirementStorage> requirements;
    core::PlannedBuffer<char> request_text;
    core::PlannedBuffer<SiEnclosure> converted;
    core::PlannedBuffer<std::uint64_t> node_to_row;
    SiResolvedDefaults defaults{false, false, false};
    SiMemoryReport report{};
};
struct SiContext {
    const SiLimits& limits;
    core::ResourcePlan& resources;
    SiDiagnostic* slot;
    si::Work work;
    core::PlannedBuffer<si::ProgramScratch> program;
    std::optional<SourceIdentity> source;
    std::optional<SiQuantityKey> key;
    const InputDocumentData* snapshot = nullptr;
    std::span<const SchemaNode> nodes;
    std::optional<std::uint64_t> current;
    std::optional<SchemaDimension> dimensions;
    std::optional<SiEnclosure> achieved;
    std::optional<SiRequirementView> requirement;
    std::uint64_t owned = 0, scratch = 0, start_charged;
    std::uint64_t automatic_scratch = 0, exact_scratch_peak = 0, opaque_exact_work = 0;
    std::uint64_t planned_owned = 0, planned_scratch = 0, planned_work = 0, admission = 0;
    SiContext(const SiLimits&, core::ResourcePlan&, SiDiagnostic*);
    [[nodiscard]] std::uint64_t limit(SiQuota) const;
    [[nodiscard]] std::uint64_t add(std::uint64_t, std::uint64_t);
    [[nodiscard]] std::uint64_t multiply(std::uint64_t, std::uint64_t);
    void require(SiQuota, std::uint64_t);
    void admit(std::uint64_t, std::uint64_t, std::uint64_t, std::uint64_t,
               std::optional<std::uint64_t> allocated_scratch = {},
               std::uint64_t cumulative_extra = 0, std::uint64_t diagnostic_request_bytes = 0);
    void finish();
    [[noreturn]] void reject(SiReason, core::FailureCode = core::FailureCode::InvalidInput,
                             std::optional<SiQuota> = {}, std::optional<std::uint64_t> = {});
    template <class T>
    [[nodiscard]] core::PlannedBuffer<T> allocate(std::uint64_t n, bool temporary = false) {
        const auto bytes = multiply(n, sizeof(T));
        if (n > static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max()))
            reject(SiReason::Overflow, core::FailureCode::ResourceExceeded);
        if (temporary) {
            scratch = add(scratch, bytes);
            require(SiQuota::ScratchBytes, scratch);
        } else {
            owned = add(owned, bytes);
            require(SiQuota::OwnedBytes, owned);
        }
        work.step(bytes);
        auto result =
            resources.allocate<T>(static_cast<core::Index>(n), core::Operation::PrepareModel);
        if (!result)
            throw result.error();
        require(SiQuota::CumulativeBytes, resources.allocated_bytes() - start_charged);
        return std::move(*result);
    }
    [[nodiscard]] SiMemoryReport report() const noexcept {
        return {owned,
                scratch + automatic_scratch + exact_scratch_peak,
                resources.allocated_bytes() - start_charged,
                work.used,
                planned_owned,
                planned_scratch,
                planned_work,
                admission};
    }
};
struct SiFactory {
    [[nodiscard]] static bool owns_key(const SiCatalogData& owner,
                                       const SiQuantityKey& key) noexcept {
        return owner.key_owner == key.owner_;
    }
    [[nodiscard]] static const SchemaData* storage(const SchemaCase& value) noexcept {
        return value.data_.size() ? &value.data_.values()[0] : nullptr;
    }
    [[nodiscard]] static const SiCatalogData* storage(const SiCatalog& value) noexcept {
        return value.data_.size() ? &value.data_.values()[0] : nullptr;
    }
    [[nodiscard]] static const SiRequestData* storage(const SiRequest& value) noexcept {
        return value.data_.size() ? &value.data_.values()[0] : nullptr;
    }
    [[nodiscard]] static const SiData* storage(const SiCase& value) noexcept {
        return value.data_.size() ? &value.data_.values()[0] : nullptr;
    }
    [[nodiscard]] static SiQuantityKey key(std::shared_ptr<const SiKeyOwner> owner,
                                           SchemaIdentity identity, std::uint64_t node,
                                           std::uint16_t type, std::uint16_t role) noexcept {
        return {std::move(owner), identity, node, type, role};
    }
    [[nodiscard]] static std::uint64_t node(const SiQuantityKey& key) noexcept {
        return key.node_;
    }
    [[nodiscard]] static SiCatalog catalog(core::PlannedBuffer<SiCatalogData> p) noexcept {
        return SiCatalog{std::move(p)};
    }
    [[nodiscard]] static SiRequest request(core::PlannedBuffer<SiRequestData> p) noexcept {
        return SiRequest{std::move(p)};
    }
    [[nodiscard]] static SiCase result(core::PlannedBuffer<SiData> p) noexcept {
        return SiCase{std::move(p)};
    }
    static void diagnostic(SiDiagnostic& slot, const SiDiagnosticDetail& fixed,
                           core::PlannedBuffer<char> path = {}) noexcept {
        slot.detail_ = fixed;
        slot.path_ = std::move(path);
    }
    static void diagnostic_requirement(SiDiagnostic& slot, core::PlannedBuffer<char> text,
                                       SiRequirementView view) noexcept {
        slot.requirement_text_ = std::move(text);
        slot.requirement_ = view;
    }
};
[[nodiscard]] bool same_identity(const SchemaIdentity&, const SchemaIdentity&) noexcept;
[[nodiscard]] std::uint64_t si_snapshot_bytes(SiContext&, const InputDocumentData&);
void si_clone_snapshot(SiContext&, const InputDocumentData&, InputDocumentData&);
[[nodiscard]] std::uint64_t si_schema_clone_bytes(SiContext&, const SchemaData&);
[[nodiscard]] SchemaCase si_clone_schema(SiContext&, const SchemaData&);
[[nodiscard]] SiUnits si_units(const SchemaData&);
[[nodiscard]] bool si_numeric_node(const SchemaNode&) noexcept;
[[nodiscard]] SiEntryStorage si_entry(SiContext&, std::span<const SchemaNode>, std::uint64_t);
[[nodiscard]] SiRecipeView si_recipe(const InputDocumentData&, std::span<const SchemaNode>,
                                     std::span<const char>, std::uint64_t, SiUnits);
[[nodiscard]] SiCatalogEntry si_catalog_entry(const SiCatalogData&, const SchemaIdentity&,
                                              std::uint64_t);
[[nodiscard]] si::Interval si_exact_interval(SiContext&, const SiExactInput&);
[[nodiscard]] si::Interval si_unit_interval(SiContext&, const SiRecipeView&,
                                            const si::Interval& pi);
[[nodiscard]] std::optional<SiEnclosure> si_exact_zero_center(SiContext&, const SiRecipeView&,
                                                              bool);
[[nodiscard]] SiRequirementView si_requirement_view(const SiRequirementStorage&,
                                                    std::span<const char>);
[[nodiscard]] SiEnclosure si_evaluate_quantity(SiContext&, const SiRecipeView&,
                                               const SiRequirementView&, const si::Interval&);
[[nodiscard]] core::Result<std::uint64_t> si_decimal_token_size(SiDecimalView) noexcept;
} // namespace uvlas::io::detail
