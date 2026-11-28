#include "observation_file.hpp"
#include "source_digest.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <fstream>
#include <iterator>
#include <string>

namespace {
using namespace uvlas;
using namespace uvlas::io;
using namespace uvlas::io::detail;
using structure::ExactRational;

class Reader final : public ReferenceFileReader {
public:
    explicit Reader(std::span<const std::byte> bytes) : bytes_(bytes), length_(bytes.size()) {}

    core::Result<ReferenceFileInfo> info() noexcept override {
        return ReferenceFileInfo{{1, {}}, length_, {}, true};
    }

    core::Result<std::size_t> read(std::span<std::byte> target) noexcept override {
        const auto count = std::min(target.size(), bytes_.size());
        std::ranges::copy(bytes_.first(count), target.begin());
        bytes_ = bytes_.subspan(count);
        return count;
    }

private:
    std::span<const std::byte> bytes_;
    std::size_t length_;
};

const structure::ExactLimits exact_limits{4096,
                                          1,
                                          1,
                                          0,
                                          0,
                                          16'384,
                                          65'536,
                                          structure::ExactEstimateAdmission{1ULL << 28},
                                          structure::ExactTimingMode::Disabled};
const InputLimits json_limits{
    1U << 20, 120,  100'000,  100'000,
    100'000,  4096, 1U << 20, {ParserAllocationMode::EstimateAdmission, 1ULL << 27}};
constexpr std::array<std::uint64_t, 2> observations{1, 2};

std::uint64_t external_bytes(const ObservationFileLimits& limits) {
    return std::get<structure::ExactEstimateAdmission>(limits.exact.external_memory).bytes;
}

} // namespace

TEST_CASE("Static observation file bounds aggregate work and retains exact owned values",
          "[FD-09][V-USE][verification][negative]") {
    std::ifstream file{UVLAS_STATIC_OBSERVATION_FILE, std::ios::binary};
    REQUIRE(file);
    const std::string data{std::istreambuf_iterator<char>{file}, {}};
    Reader reader{std::as_bytes(std::span{data})};
    const auto digest =
        decode_source_digest("fb599adc52513442ab45129c0e90aadcf0bf7647fe216117d51d4563e9bbd5f4");
    const auto condition =
        decode_source_digest("378fea09ac121f7e2240224bfacb63c983ab73b695fe39199ab710685c900b0e");
    if (!digest || !condition) {
        FAIL("Expected the published SHA-256 values");
        return;
    }

    const SchemaIdentity parent{{SourceId{1}, 12, {}}, 2, 2};
    const RootSnapshotLocation location{parent, SourceId{6}, ReferenceRole::ObservationReference,
                                        42};
    core::ResourcePlan plan{1ULL << 27};
    auto snapshot = read_reference_file(reader, *digest, {1U << 20, 1U << 20, 100, 1U << 20}, plan,
                                        {}, nullptr, {});
    if (!snapshot) {
        FAIL("Could not authenticate the published reference file");
        return;
    }
    const AcquiredRootSnapshot source{location, std::move(*snapshot), {}};
    const ObservationFileLimits initial{json_limits, exact_limits, 10'000'000};
    auto available = initial;
    std::optional<ObservationFileDiagnostic> diagnostic;
    auto parsed = read_observation_file(source, parent, *condition, observations, available, plan,
                                        {}, diagnostic);
    if (!parsed) {
        FAIL("Could not parse the published reference file");
        return;
    }
    const auto work = initial.maximum_work - available.maximum_work;
    const auto bytes = external_bytes(initial) - external_bytes(available);
    REQUIRE(work > 0);
    REQUIRE(bytes > 0);
    CHECK(parsed->source.parent.source == parent.source);
    CHECK(parsed->source.parent.contract_revision == parent.contract_revision);
    CHECK(parsed->source.parent.representation_revision == parent.representation_revision);
    CHECK(parsed->source.source == location.source);
    CHECK(parsed->source_sha256 == *digest);
    CHECK(parsed->condition_sha256 == *condition);
    REQUIRE(parsed->records.size() == 2);
    const auto& point = parsed->records.values()[0];
    const auto& reaction = parsed->records.values()[1];
    const auto& point_y = point.components[1];
    const auto& force_y = reaction.components[1];
    const auto& moment_x = reaction.components[3];
    const auto& moment_z = reaction.components[5];
    if (!point_y || !force_y || !moment_x || !moment_z) {
        FAIL("Missing published component");
        return;
    }
    const auto displacement = ExactRational::from_integer_pair("1", "512", exact_limits, plan);
    const auto force =
        ExactRational::from_integer_pair("-4167633984375", "4194304", exact_limits, plan);
    const auto moment =
        ExactRational::from_integer_pair("4167633984375", "4194304", exact_limits, plan);
    if (!displacement || !force || !moment) {
        FAIL("Could not construct the analytical reference");
        return;
    }
    CHECK(point_y->value == *displacement);
    CHECK(force_y->value == *force);
    CHECK(moment_x->value == *moment);
    CHECK(moment_z->value == *moment);

    for (const bool reserve_memory : {false, true}) {
        auto tight = initial;
        if (reserve_memory)
            tight.exact.external_memory = structure::ExactEstimateAdmission{2 * bytes - 1};
        else
            tight.maximum_work = 2 * work - 1;
        const auto first = read_observation_file(source, parent, *condition, observations, tight,
                                                 plan, {}, diagnostic);
        if (!first) {
            FAIL("Expected the first file to fit within the remaining budget");
            return;
        }
        const auto second = read_observation_file(source, parent, *condition, observations, tight,
                                                  plan, {}, diagnostic);
        REQUIRE_FALSE(second);
        CHECK(second.error().code == core::FailureCode::ResourceExceeded);
        REQUIRE(first->records.size() == 2);
        const auto& retained = first->records.values()[0].components[1];
        if (!retained) {
            FAIL("The previously owned component is missing");
            return;
        }
        CHECK(retained->value == *displacement);
    }
    auto different_parent = parent;
    different_parent.source.digest.bytes[0] = std::byte{1};
    auto fresh = initial;
    auto rejected = read_observation_file(source, different_parent, *condition, observations, fresh,
                                          plan, {}, diagnostic);
    REQUIRE_FALSE(rejected);
    CHECK(rejected.error().code == core::FailureCode::StaleGeneration);
    std::stop_source cancellation;
    cancellation.request_stop();
    rejected = read_observation_file(source, parent, *condition, observations, fresh, plan,
                                     cancellation.get_token(), diagnostic);
    REQUIRE_FALSE(rejected);
    CHECK(rejected.error().code == core::FailureCode::Cancelled);
    CHECK(fresh.maximum_work == initial.maximum_work);
    CHECK(external_bytes(fresh) == external_bytes(initial));
    CHECK(read_observation_file(source, parent, *condition, observations, fresh, plan, {},
                                diagnostic));
}
