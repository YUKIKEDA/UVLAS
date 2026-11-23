#pragma once

#include <span>
#include <stop_token>
#include <string_view>
#include <uvlas/core/resource_plan.hpp>
#include <uvlas/io/si_case.hpp>
#include <uvlas/numerics/precision.hpp>

namespace uvlas::io::detail {

enum class ObservationSide : std::uint8_t { None, Interior, Lower, Upper };

/// An already-qualified component, in the order supplied to the observation evaluator.
struct ObservationCsvColumn {
    core::ObservationId observation;
    std::uint64_t component;
    std::string_view quantity;
    std::string_view unit;
    std::string_view basis;
    std::string_view reference_point;
    std::optional<core::LayerId> layer;
    ObservationSide side;
};

/// An accepted sample, with one qualified enclosure per column in the same order.
struct ObservationCsvSample {
    numerics::TwoWordValue time_seconds{};
    std::uint64_t checkpoint = 0;
    std::uint64_t segment = 0;
    std::span<const SiEnclosure> values;
};

struct ObservationCsvLimits {
    std::size_t maximum_columns;
    std::size_t maximum_samples;
    std::size_t maximum_label_bytes;
    std::size_t maximum_output_bytes;
};

/// Encode UTF-8 records into owned bytes. Views must remain stable during the call.
/// Empty series produce a header. Authority for acceptance and physical qualification
/// remains with the caller. Any validation, budget or cancellation failure returns no bytes.
[[nodiscard]] core::Result<core::PlannedBuffer<char>>
encode_observation_csv(std::span<const ObservationCsvColumn> columns,
                       std::span<const ObservationCsvSample> samples, ObservationCsvLimits limits,
                       core::ResourcePlan& resources, const std::stop_token& stop = {}) noexcept;

} // namespace uvlas::io::detail
