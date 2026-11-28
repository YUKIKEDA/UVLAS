#pragma once

#include <stop_token>
#include <uvlas/io/input_document.hpp>

namespace uvlas::io::detail {

struct ObservationConditionLimits {
    std::uint64_t maximum_work;
    std::uint64_t maximum_hash_bytes;
};

[[nodiscard]] core::Result<SourceDigest>
observation_condition_digest(const InputNodeView&, ObservationConditionLimits, core::ResourcePlan&,
                             const std::stop_token&) noexcept;

} // namespace uvlas::io::detail
