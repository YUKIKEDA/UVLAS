#include "reference_declarations_fixture.hpp"
#include "reference_input.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <chrono>
#include <fstream>
#include <type_traits>

namespace {
namespace io = uvlas::io;
namespace core = uvlas::core;
namespace rd = uvlas::io::detail;
namespace fixture = uvlas::test::declaration_fixture;
namespace fs = std::filesystem;

constexpr io::SourceId source{19};
constexpr rd::ManifestLimits manifest_limits{
    128, {fixture::budget, fixture::budget, fixture::budget, fixture::budget}};
constexpr std::string_view name = "\xc3\xa9.json";
constexpr char8_t native_name[] = u8"\u00e9.json";

static_assert(!std::is_constructible_v<rd::ReferenceInput, io::SchemaCase, rd::ReadRoot,
                                       rd::ReferenceFileId>);
static_assert(!std::is_copy_constructible_v<rd::ReferenceInput>);

struct Tree {
    fs::path root;

    Tree() {
        const auto tick = std::chrono::steady_clock::now().time_since_epoch().count();
        for (int attempt = 0; attempt < 16; ++attempt) {
            root = fs::temp_directory_path() / ("uvlas-reference-input-" + std::to_string(tick) +
                                                "-" + std::to_string(attempt));
            if (fs::create_directory(root)) {
                return;
            }
        }
        FAIL("Could not create a private test directory");
    }
    ~Tree() {
        std::error_code ignore;
        fs::remove_all(root, ignore);
    }
    Tree(const Tree&) = delete;
    Tree& operator=(const Tree&) = delete;
};

void write_file(const fs::path& target, std::string_view bytes) {
    std::ofstream out(target, std::ios::binary);
    out << bytes;
    out.close();
    REQUIRE(out.good());
}

std::span<const std::byte> bytes_of(std::string_view text) noexcept {
    return std::as_bytes(std::span{text});
}

rd::ReferenceInputLimits limits() {
    const io::InputLimits input{
        10'000'000, 256,
        1'000'000,  1'000'000,
        10'000'000, 10'000'000,
        10'000'000, {io::ParserAllocationMode::EstimateAdmission, fixture::budget}};
    std::array<std::uint64_t, static_cast<std::size_t>(io::SchemaQuota::Count)> caps{};
    caps.fill(fixture::budget);
    return {input, io::SchemaLimits{caps}};
}

constexpr std::string_view static_case() {
    for (const auto& example : uvlas::test::case_schema_corpus::roots) {
        if (example.name == "structural_static") {
            return example.json;
        }
    }
    return {};
}

consteval io::SourceDigest digest_literal(std::string_view hex) {
    constexpr std::string_view digits = "0123456789abcdef";
    io::SourceDigest digest{};
    for (std::size_t i = 0; i < digest.bytes.size(); ++i) {
        digest.bytes[i] =
            static_cast<std::byte>(16 * digits.find(hex[2 * i]) + digits.find(hex[2 * i + 1]));
    }
    return digest;
}

// Python 3.12 hashlib SHA-256 over both concatenated structural_static corpus literals.
constexpr auto static_case_digest =
    digest_literal("ecfba5a3cdf129b195b7e1504c74c2981852f6ee57617ac320ea14420ffd9a38");
constexpr auto abc_digest =
    digest_literal("ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");

void check_failure(const core::Result<rd::ReferenceInput>& result, core::FailureCode code) {
    REQUIRE_FALSE(result);
    CHECK(result.error().operation == core::Operation::PrepareModel);
    CHECK(result.error().code == code);
}

#ifdef __linux__
std::size_t held_descriptors(const fs::path& path) {
    std::size_t result = 0;
    for (const auto& entry : fs::directory_iterator("/proc/self/fd")) {
        std::error_code ignore;
        if (fs::read_symlink(entry.path(), ignore) == path) {
            ++result;
        }
    }
    return result;
}
#endif

} // namespace

TEST_CASE("File and Memory reference inputs retain one checked case and their actual origins",
          "[IO4-01][N4-02][integration][property]") {
    static_assert(static_case().size() == 8611);
    const Tree tree;
    const auto sub = tree.root / "child";
    REQUIRE(fs::create_directory(sub));
    write_file(tree.root / native_name, "not the manifest");
    write_file(sub / native_name, static_case());
    write_file(tree.root / "copy.bin", "abd");
    write_file(sub / "copy.bin", "abc");
    fs::create_hard_link(sub / native_name, sub / "another-name.json");

    {
        auto file = [&] {
            core::ResourcePlan plan{fixture::budget};
            rd::ReferenceInputDiagnostic diagnostic;
            auto result = rd::ReferenceInput::from_file(sub, name, source, manifest_limits,
                                                        limits(), plan, {}, &diagnostic);
            CHECK_FALSE(diagnostic.stage);
            return fixture::take(std::move(result));
        }();
        auto owned = std::move(file);
        // NOLINTNEXTLINE(bugprone-use-after-move): moved-from rejection is part of the contract.
        CHECK_FALSE(file.view());

        const auto view = fixture::take(owned.view());
        REQUIRE(view.root != nullptr);
        REQUIRE(view.manifest_identity);
        const auto identity = fixture::take(view.schema.identity());
        CHECK(identity.source.source == source);
        CHECK(identity.source.byte_count == static_case().size());
        CHECK(identity.source.digest == static_case_digest);

        auto memory = [&] {
            std::string buffer{static_case()};
            core::ResourcePlan plan{fixture::budget};
            auto result = fixture::take(rd::ReferenceInput::from_memory(
                bytes_of(buffer), source, std::nullopt, limits(), plan, {}, nullptr));
            buffer.assign(buffer.size(), '!');
            return result;
        }();
        const auto memory_view = fixture::take(memory.view());
        CHECK(memory_view.root == nullptr);
        CHECK_FALSE(memory_view.manifest_identity);
        CHECK(fixture::same_identity(identity, fixture::take(memory_view.schema.identity())));

#ifdef __linux__
        CHECK(held_descriptors(sub) == 1);
        CHECK(held_descriptors(sub / native_name) == 1);
#else
        std::error_code manifest_rename, root_rename;
        fs::rename(sub / native_name, sub / "renamed.json", manifest_rename);
        fs::rename(sub, tree.root / "renamed-child", root_rename);
        REQUIRE(manifest_rename);
        REQUIRE(root_rename);
#endif

        core::ResourcePlan plan{fixture::budget};
        const auto alias_path = fixture::take(
            rd::normalize_reference_path("another-name.json", {64, 64, 1}, plan, {}, nullptr));
        auto alias = fixture::take(view.root->open_file(alias_path, plan, {}));
        CHECK(view.manifest_identity == fixture::take(alias.info()).identity);

        const auto copy_path =
            fixture::take(rd::normalize_reference_path("copy.bin", {64, 64, 1}, plan, {}, nullptr));
        auto copy = fixture::take(view.root->open_file(copy_path, plan, {}));
        const auto snapshot = fixture::take(
            rd::read_reference_file(copy, abc_digest, {3, 4, 2, 3}, plan, {}, nullptr));
        if (const auto payload = snapshot.view()) {
            CHECK(payload->digest == abc_digest);
            CHECK(std::ranges::equal(payload->bytes, bytes_of("abc")));
        } else {
            FAIL("Read did not retain the external payload");
        }
    }
#ifdef __linux__
    CHECK(held_descriptors(sub) == 0);
    CHECK(held_descriptors(sub / native_name) == 0);
#else
    fs::rename(sub / native_name, sub / "renamed.json");
    fs::rename(sub, tree.root / "renamed-child");
#endif
}

TEST_CASE("Memory reference input requires a root for the eight declared roles",
          "[IO4-01][N4-02][integration][negative]") {
    const Tree tree;
    const auto supplied = tree.root / "supplied";
    REQUIRE(fs::create_directory(supplied));
    write_file(supplied / "copy.bin", "abc");
    const auto wire = fixture::all_roles().wire();
    const auto bytes = bytes_of(wire);
    core::ResourcePlan plan{fixture::budget};
    rd::ReferenceInputDiagnostic diagnostic;
    check_failure(rd::ReferenceInput::from_memory(bytes, source, std::nullopt, limits(), plan, {},
                                                  &diagnostic),
                  core::FailureCode::InvalidInput);
    CHECK(diagnostic.stage == rd::ReferenceInputStage::Root);
    check_failure(
        rd::ReferenceInput::from_memory(bytes, source, std::nullopt, limits(), plan, {}, nullptr),
        core::FailureCode::InvalidInput);

    auto held = fixture::take(rd::ReadRoot::open(supplied, {}));
#ifdef __linux__
    const auto relocated = tree.root / "relocated";
    fs::rename(supplied, relocated);
    REQUIRE(fs::create_directory(supplied));
    write_file(supplied / "copy.bin", "abd");
#endif
    const auto input = fixture::take(rd::ReferenceInput::from_memory(
        bytes, source, std::move(held), limits(), plan, {}, &diagnostic));
    CHECK_FALSE(diagnostic.stage);
    const auto view = fixture::take(input.view());
    REQUIRE(view.root != nullptr);
    CHECK_FALSE(view.manifest_identity);
    const auto copy_path =
        fixture::take(rd::normalize_reference_path("copy.bin", {64, 64, 1}, plan, {}, nullptr));
    auto copy = fixture::take(view.root->open_file(copy_path, plan, {}));
    const auto snapshot =
        fixture::take(rd::read_reference_file(copy, abc_digest, {3, 4, 2, 3}, plan, {}, nullptr));
    if (const auto payload = snapshot.view()) {
        CHECK(std::ranges::equal(payload->bytes, bytes_of("abc")));
    } else {
        FAIL("Read did not retain the external payload");
    }

    const auto declarations = fixture::take(
        rd::collect_case_file_references(view.schema, fixture::limits(), plan, {}, nullptr));
    CHECK(fixture::take(declarations.size()) == 8);
    CHECK(fixture::take(declarations.identity()).source.source == source);

    // NOLINTNEXTLINE(bugprone-use-after-move): a moved-from root cannot qualify an input.
    auto closed = std::optional{std::move(held)};
    check_failure(rd::ReferenceInput::from_memory(bytes, source, std::move(closed), limits(), plan,
                                                  {}, &diagnostic),
                  core::FailureCode::InvalidInput);
    CHECK(diagnostic.stage == rd::ReferenceInputStage::Root);
}

TEST_CASE("Reference acquisition rejects invalid names, bytes, quota and cancellation",
          "[IO4-01][N4-08][integration][negative]") {
    const Tree tree;
    core::ResourcePlan plan{fixture::budget};
    rd::ReferenceInputDiagnostic diagnostic;
    for (const std::string_view bad :
         {"", ".", "..", "dir/file", "dir\\file", "/absolute", "a:b", "p/../file", "\xc0\x80",
          "\xe0\x9f\xbf", "\xed\xa0\x80", "\xf4\x90\x80\x80"}) {
        CAPTURE(bad);
        check_failure(rd::ReferenceInput::from_file(tree.root, bad, source, manifest_limits,
                                                    limits(), plan, {}, &diagnostic),
                      core::FailureCode::InvalidInput);
        CHECK(diagnostic.stage == rd::ReferenceInputStage::ManifestName);
    }
    check_failure(rd::ReferenceInput::from_file(tree.root, name, source, manifest_limits, limits(),
                                                plan, {}, &diagnostic),
                  core::FailureCode::InvalidInput);
    CHECK(diagnostic.stage == rd::ReferenceInputStage::ManifestOpen);

    write_file(tree.root / native_name, "{");
    check_failure(rd::ReferenceInput::from_file(tree.root, name, source, manifest_limits, limits(),
                                                plan, {}, &diagnostic),
                  core::FailureCode::InvalidInput);
    CHECK(diagnostic.stage == rd::ReferenceInputStage::Document);
    CHECK(diagnostic.document.has_value());

    auto broken = fixture::root();
    broken.field("model").erase("materials");
    const auto wire = broken.wire();
    check_failure(rd::ReferenceInput::from_memory(bytes_of(wire), source, std::nullopt, limits(),
                                                  plan, {}, &diagnostic),
                  core::FailureCode::InvalidInput);
    CHECK(diagnostic.stage == rd::ReferenceInputStage::Schema);
    CHECK(diagnostic.schema.detail().has_value());
    check_failure(rd::ReferenceInput::from_memory(bytes_of(wire), source, std::nullopt, limits(),
                                                  plan, {}, nullptr),
                  core::FailureCode::InvalidInput);

    write_file(tree.root / native_name, static_case());
    auto insufficient = manifest_limits;
    insufficient.file.maximum_file_bytes = static_case().size() - 1;
    check_failure(rd::ReferenceInput::from_file(tree.root, name, source, insufficient, limits(),
                                                plan, {}, &diagnostic),
                  core::FailureCode::ResourceExceeded);
    CHECK(diagnostic.stage == rd::ReferenceInputStage::ManifestRead);
    if (const auto file = diagnostic.file) {
        CHECK(file->reason == rd::ReferenceFileFailure::FileBytes);
    } else {
        FAIL("File quota rejection did not retain its reason");
    }

    // NOLINTNEXTLINE(misc-const-correctness): MSVC 14.51 request_stop() is non-const.
    std::stop_source stopped;
    REQUIRE(stopped.request_stop());
    const auto before = plan.allocated_bytes();
    check_failure(rd::ReferenceInput::from_file(tree.root, name, source, manifest_limits, limits(),
                                                plan, stopped.get_token(), &diagnostic),
                  core::FailureCode::Cancelled);
    CHECK_FALSE(diagnostic.stage);
    check_failure(rd::ReferenceInput::from_memory(bytes_of(static_case()), source, std::nullopt,
                                                  limits(), plan, stopped.get_token(), &diagnostic),
                  core::FailureCode::Cancelled);
    CHECK_FALSE(diagnostic.stage);
    CHECK(plan.allocated_bytes() == before);
}

TEST_CASE("Both acquisition routes account for their complete allocation on one plan",
          "[IO4-01][N4-02][integration][negative]") {
    const Tree tree;
    write_file(tree.root / native_name, static_case());
    const auto check_budget = [](const auto& acquire) {
        core::ResourcePlan probe{fixture::budget};
        const auto input = fixture::take(acquire(probe));
        const auto view = fixture::take(input.view());
        const auto schema = fixture::take(view.schema.memory_report());
        REQUIRE(probe.allocated_bytes() >= schema.cumulative_bytes);
        REQUIRE(schema.admission_bytes >= schema.cumulative_bytes);
        // Schema admits a failure reserve as well as charging successful allocations.
        // The acquisition and parsing charges must still be present at schema admission.
        const auto required =
            probe.allocated_bytes() - schema.cumulative_bytes + schema.admission_bytes;
        REQUIRE(required > 0);

        core::ResourcePlan exact{required};
        REQUIRE(acquire(exact));
        CHECK(exact.allocated_bytes() == probe.allocated_bytes());

        core::ResourcePlan one_short{required - 1};
        check_failure(acquire(one_short), core::FailureCode::ResourceExceeded);
    };
    check_budget([&](core::ResourcePlan& plan) {
        return rd::ReferenceInput::from_file(tree.root, name, source, manifest_limits, limits(),
                                             plan, {}, nullptr);
    });
    check_budget([](core::ResourcePlan& plan) {
        return rd::ReferenceInput::from_memory(bytes_of(static_case()), source, std::nullopt,
                                               limits(), plan, {}, nullptr);
    });
}
