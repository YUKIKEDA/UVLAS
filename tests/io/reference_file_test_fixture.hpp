#pragma once

#include "reference_declarations_fixture.hpp"
#include "reference_input.hpp"
#include "reference_root_plan.hpp"

#include <catch2/catch_test_macros.hpp>

#include <chrono>
#include <fstream>

namespace uvlas::test::reference_file_fixture {
namespace fixture = uvlas::test::declaration_fixture;
namespace rd = uvlas::io::detail;
namespace fs = std::filesystem;

constexpr rd::ManifestLimits manifest_limits{
    128, {fixture::budget, fixture::budget, fixture::budget, fixture::budget}};
constexpr rd::RootReferenceLimits root_limits{
    fixture::budget, fixture::budget, fixture::budget, fixture::budget, 64, fixture::budget, 0};

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

inline void write_file(const fs::path& target, std::string_view bytes) {
    std::ofstream out(target, std::ios::binary);
    out << bytes;
    out.close();
    REQUIRE(out.good());
}

inline std::span<const std::byte> bytes_of(std::string_view text) noexcept {
    return std::as_bytes(std::span{text});
}

inline rd::ReferenceInputLimits limits() {
    const io::InputLimits input{
        10'000'000, 256,
        1'000'000,  1'000'000,
        10'000'000, 10'000'000,
        10'000'000, {io::ParserAllocationMode::EstimateAdmission, fixture::budget}};
    std::array<std::uint64_t, static_cast<std::size_t>(io::SchemaQuota::Count)> caps{};
    caps.fill(fixture::budget);
    return {input, io::SchemaLimits{caps}};
}

inline void declare_path(fixture::JsonToken& file, std::string_view path, std::string_view hash) {
    file.set("path", fixture::string(path));
    file.set("sha256", fixture::string(hash));
}

inline rd::PreparedRootReferences prepare(const rd::ReferenceInput& input) {
    core::ResourcePlan plan{fixture::budget};
    const auto view = fixture::take(input.view());
    auto declared = fixture::take(
        rd::collect_case_file_references(view.schema, fixture::limits(), plan, {}, nullptr));
    return fixture::take(
        rd::prepare_root_references(std::move(declared), root_limits, plan, {}, nullptr));
}

} // namespace uvlas::test::reference_file_fixture
