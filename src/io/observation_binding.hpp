#pragma once

#include "observation_file.hpp"

namespace uvlas::io::detail {

struct ObservationBindingDiagnostic {
    RootSnapshotLocation declaration;
    std::optional<core::ObservationId> observation;
    std::optional<core::Index> case_index;
    std::optional<InputNodeOrdinal> input_node;
};

[[nodiscard]] core::Result<void>
bind_observation_file(const SchemaCase&, const ObservationFile&, const structure::ExactRational&,
                      structure::ExactLimits&, std::uint64_t&, core::ResourcePlan&,
                      const std::stop_token&,
                      std::optional<ObservationBindingDiagnostic>&) noexcept;

} // namespace uvlas::io::detail
