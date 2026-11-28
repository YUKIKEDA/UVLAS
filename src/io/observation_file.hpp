#pragma once

#include "reference_root_snapshot.hpp"

#include <array>
#include <uvlas/structure/exact_polynomial.hpp>

namespace uvlas::io::detail {

struct ObservationFileLimits {
    ObservationFileLimits(InputLimits syntax, structure::ExactLimits arithmetic,
                          std::uint64_t work) noexcept
        : json(syntax), exact(arithmetic), maximum_work(work) {}

    InputLimits json;
    structure::ExactLimits exact;
    std::uint64_t maximum_work;
};

struct ObservationFileScalar {
    structure::ExactRational value, uncertainty;
};

enum class ObservationPointSide : std::uint8_t { Interior, Lower, Upper };
struct ObservationFilePoint {
    std::uint64_t shell;
    ObservationPointSide side;
    std::array<structure::ExactRational, 3> coordinate;
};

struct ObservationFileReaction {
    std::uint64_t support;
    std::array<structure::ExactRational, 3> origin_metres;
};

struct ObservationFileRecord {
    core::ObservationId id{0};
    std::optional<InputNodeOrdinal> input_node;
    std::variant<std::monostate, ObservationFilePoint, ObservationFileReaction> descriptor;
    structure::ExactRational time_seconds;
    std::array<std::optional<ObservationFileScalar>, 6> components;
};

struct ObservationFileDiagnostic {
    RootSnapshotLocation declaration;
    std::optional<InputNodeOrdinal> node;
    std::optional<std::size_t> byte_offset;
};

struct ObservationFile {
    RootSnapshotLocation source;
    SourceDigest source_sha256;
    SourceDigest condition_sha256;
    core::PlannedBuffer<ObservationFileRecord> records;
};

[[nodiscard]] core::Result<ObservationFile>
read_observation_file(const AcquiredRootSnapshot&, const SchemaIdentity&, const SourceDigest&,
                      std::span<const std::uint64_t>, ObservationFileLimits&, core::ResourcePlan&,
                      const std::stop_token&, std::optional<ObservationFileDiagnostic>&) noexcept;

} // namespace uvlas::io::detail
