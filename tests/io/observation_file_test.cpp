#include "observation_file.hpp"
#include "reference_file_test_fixture.hpp"
#include "reference_resolution.hpp"
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

#if defined(__linux__) || defined(_WIN32)
TEST_CASE("Observation graph publishes only after both sources succeed",
          "[FD-09][V-USE][verification][negative]") {
    namespace files = uvlas::test::reference_file_fixture;
    namespace fixture = uvlas::test::declaration_fixture;

    const auto read = [](const char* path) {
        std::ifstream file{path, std::ios::binary};
        REQUIRE(file);
        return fixture::json(std::string{std::istreambuf_iterator<char>{file}, {}});
    };
    auto manifest = read(UVLAS_STATIC_OBSERVATION_CASE);
    const auto shared = read(UVLAS_STATIC_OBSERVATION_FILE);
    auto& declarations = manifest.field("observations").elements;
    const auto& records = shared.field("observations").elements;
    REQUIRE(declarations.size() == 2);
    REQUIRE(records.size() == 2);
    const files::Tree tree;
    const std::array names{"point.json", "reaction.json"};
    std::array<std::string, 2> payload;
    for (std::size_t i = 0; i < names.size(); ++i) {
        auto separated = shared;
        separated.field("observations").elements = {records[i]};
        payload[i] = separated.wire();
        const auto digest = source_digest(files::bytes_of(payload[i]));
        constexpr std::string_view digits = "0123456789abcdef";
        std::string hex;
        for (const auto byte : digest.bytes) {
            const auto x = std::to_integer<unsigned>(byte);
            hex.push_back(digits[x >> 4U]);
            hex.push_back(digits[x & 15U]);
        }
        auto& source = declarations[i].field("comparison").field("reference").field("source");
        source.set("source_id", fixture::string(std::to_string(10 + i)));
        files::declare_path(source, names[i], hex);
        files::write_file(tree.root / names[i], payload[i]);
    }
    files::write_file(tree.root / "case.json", manifest.wire());
    core::ResourcePlan input_plan{fixture::budget};
    const auto input = fixture::take(
        ReferenceInput::from_file(tree.root, "case.json", SourceId{1}, files::manifest_limits,
                                  files::limits(), input_plan, {}, nullptr));

    constexpr RootInventoryLimits inventory{
        {{1U << 20, 1U << 20, 128, 1U << 20}, 8, 128}, 2, 2U << 20, 2U << 20, 256, 2U << 20, 8};
    constexpr DependencyWalkLimits walk{2, 2, 1, 10'000, 128};
    ObservationLeafResults observations;
    RootSnapshotDiagnostic acquisition;
    LeafReferenceDiagnostic diagnostic;
    const auto resolve = [&](const std::stop_token& stop = {}) {
        auto prepared = files::prepare(input);
        core::ResourcePlan plan{fixture::budget};
        return resolve_leaf_references(
            input, std::move(prepared), inventory, walk,
            {{},
             {},
             {},
             ObservationLeafContext{
                 {1U << 20, 1U << 20}, {json_limits, exact_limits, 10'000'000}, 2, &observations}},
            plan, stop, &acquisition, &diagnostic);
    };

    files::write_file(tree.root / names[1], payload[1] + " ");
    const auto rejected = resolve();
    REQUIRE_FALSE(rejected);
    REQUIRE(acquisition.location);
    CHECK(acquisition.location->source == SourceId{11});
    CHECK(observations.populated == 0);
    CHECK(observations.files.size() == 0);

    files::write_file(tree.root / names[1], payload[1]);
    REQUIRE(resolve());
    REQUIRE(observations.populated == 2);
    const auto* const published = observations.files.values().data();
    REQUIRE(published[0]);
    REQUIRE(published[1]);
    CHECK(published[0]->source.source == SourceId{10});
    CHECK(published[1]->source.source == SourceId{11});

    files::write_file(tree.root / names[1], payload[1] + " ");
    REQUIRE_FALSE(resolve());
    CHECK(observations.populated == 2);
    CHECK(observations.files.values().data() == published);

    std::stop_source stopped;
    stopped.request_stop();
    const auto cancelled = resolve(stopped.get_token());
    REQUIRE_FALSE(cancelled);
    CHECK(cancelled.error().code == core::FailureCode::Cancelled);
    CHECK(observations.populated == 2);
    CHECK(observations.files.values().data() == published);

    files::write_file(tree.root / names[1], payload[1]);
    REQUIRE(resolve());
    REQUIRE(observations.populated == 2);
    CHECK(observations.files.values()[0]->source.source == SourceId{10});
    CHECK(observations.files.values()[1]->source.source == SourceId{11});
}
#endif
