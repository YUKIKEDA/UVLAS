#pragma once

#include "nastran_shell_cards.hpp"

namespace uvlas::io::detail {

struct NastranRigidGridReference {
    NastranGridId id{0};
    NastranBulkLexeme field;
    std::size_t ordinal = 0;
};

struct NastranRigidThermalField {
    NastranRealField value;
    NastranBulkLexeme field;
};

struct NastranRigidLink {
    NastranElementId id{0};
    NastranRigidGridReference independent;
    std::uint8_t dependent_components = 0;
    NastranBulkLexeme card, id_field, components_field;
    std::size_t first_dependent = 0, dependent_count = 0;
    std::optional<NastranRigidThermalField> alpha, reference_temperature;
};

enum class NastranRigidLinkFailure : std::uint8_t {
    Input,
    MissingField,
    InvalidField,
    UnsupportedThermalForm,
    UnknownGrid,
    RepeatedGrid,
    ExtraField,
    Limit
};

struct NastranRigidLinkDiagnostic {
    NastranRigidLinkFailure reason;
    NastranBulkPosition position;
};

struct NastranRigidLinkLimits {
    std::uint64_t links;
    std::uint64_t dependent_grids;
    std::uint64_t work;
};

struct NastranRigidLinkSource {
    RootSnapshotLocation location;
    SourceDigest digest;
};

struct NastranRigidLinkView {
    NastranRigidLinkSource source;
    std::span<const NastranRigidLink> links;
    std::span<const NastranRigidGridReference> dependents;
    std::uint64_t work;
};

class NastranRigidLinkCards {
public:
    NastranRigidLinkCards(NastranRigidLinkCards&&) noexcept;
    NastranRigidLinkCards& operator=(NastranRigidLinkCards&&) = delete;

    [[nodiscard]] static core::Result<NastranRigidLinkCards>
    prepare(const NastranShellCards&, NastranRigidLinkLimits, core::ResourcePlan&,
            const std::stop_token&, std::optional<NastranRigidLinkDiagnostic>* = nullptr) noexcept;

    [[nodiscard]] std::optional<NastranRigidLinkView> view() const& noexcept;
    [[nodiscard]] std::optional<NastranRigidLinkView> view() const&& = delete;

private:
    NastranRigidLinkCards(NastranRigidLinkSource, core::PlannedBuffer<NastranRigidLink>&&,
                          core::PlannedBuffer<NastranRigidGridReference>&&, std::uint64_t) noexcept;

    std::optional<NastranRigidLinkSource> source_;
    core::PlannedBuffer<NastranRigidLink> links_;
    core::PlannedBuffer<NastranRigidGridReference> dependents_;
    std::uint64_t work_;
};

} // namespace uvlas::io::detail
