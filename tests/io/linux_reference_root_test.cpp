#ifdef __linux__

#include "linux_reference_root.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <fcntl.h>
#include <filesystem>
#include <fstream>
#include <string>
#include <string_view>
#include <sys/stat.h>
#include <vector>

namespace {
namespace fs = std::filesystem;
namespace core = uvlas::core;
using namespace uvlas::io::detail;

struct Tree {
    fs::path directory;

    Tree() {
        auto pattern = (fs::temp_directory_path() / "uvlas-reference-XXXXXX").string();
        const auto* created = ::mkdtemp(pattern.data());
        REQUIRE(created != nullptr);
        directory = created;
    }
    ~Tree() {
        std::error_code ignored;
        fs::remove_all(directory, ignored);
    }
    Tree(const Tree&) = delete;
    Tree& operator=(const Tree&) = delete;
};

void write_file(const fs::path& path, std::string_view bytes) {
    std::ofstream out(path, std::ios::binary);
    out << bytes;
    out.close();
    REQUIRE(out.good());
}

[[nodiscard]] auto path_for(std::string_view input, core::ResourcePlan& plan) {
    return normalize_reference_path(input, {256, 256, 8}, plan, {}, nullptr);
}

[[nodiscard]] std::vector<int> descriptors_for(const fs::path& target) {
    std::vector<int> descriptors;
    // Procfs exposes the independent kernel descriptor view without changing the reader API.
    for (const auto& entry : fs::directory_iterator("/proc/self/fd")) {
        std::error_code ignored;
        if (fs::read_symlink(entry.path(), ignored) == target) {
            descriptors.push_back(std::stoi(entry.path().filename().string()));
        }
    }
    return descriptors;
}

void check_one_cloexec(const fs::path& target) {
    const auto descriptors = descriptors_for(target);
    REQUIRE(descriptors.size() == 1);
    const auto flags = ::fcntl(descriptors.front(), F_GETFD);
    REQUIRE(flags >= 0);
    CHECK((flags & FD_CLOEXEC) != 0);
}

consteval uvlas::io::SourceDigest digest_literal(std::string_view hex) {
    constexpr std::string_view digits = "0123456789abcdef";
    uvlas::io::SourceDigest digest{};
    for (std::size_t i = 0; i < digest.bytes.size(); ++i) {
        digest.bytes[i] =
            static_cast<std::byte>(16 * digits.find(hex[2 * i]) + digits.find(hex[2 * i + 1]));
    }
    return digest;
}

// Independent hashlib SHA-256 over the exact three input bytes.
constexpr auto abc_digest =
    digest_literal("ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
} // namespace

TEST_CASE("Linux reference lookup stays anchored after the original root path is replaced",
          "[IO4-01][verification][property]") {
    const Tree tree;
    const auto original = tree.directory / "root";
    const auto renamed = tree.directory / "renamed";
    REQUIRE(fs::create_directories(original / "sub"));
    write_file(original / "sub" / "value", "abc");
    const auto snapshot = [&] {
        const auto root = LinuxReferenceRoot::open(original, {});
        REQUIRE(root.has_value());
        check_one_cloexec(original);
        fs::rename(original, renamed);
        REQUIRE(fs::create_directories(original / "sub"));
        write_file(original / "sub" / "value", "xyz");
        const auto file_path = renamed / "sub" / "value";
        fs::create_hard_link(file_path, renamed / "alias");
        write_file(renamed / "copy", "abc");
        core::ResourcePlan plan{4096};
        const auto name = path_for("sub/./value", plan);
        const auto alias = path_for("alias", plan);
        const auto copy = path_for("copy", plan);
        REQUIRE(name.has_value());
        REQUIRE(alias.has_value());
        REQUIRE(copy.has_value());
        auto primary = root->open_file(*name, plan, {});
        REQUIRE(primary.has_value());
        check_one_cloexec(file_path);
        auto hard_link = root->open_file(*alias, plan, {});
        auto equivalent = root->open_file(*copy, plan, {});
        REQUIRE(hard_link.has_value());
        REQUIRE(equivalent.has_value());
        const auto first = primary->info();
        const auto linked = hard_link->info();
        const auto separate = equivalent->info();
        REQUIRE(first.has_value());
        REQUIRE(linked.has_value());
        REQUIRE(separate.has_value());
        CHECK(first->identity == linked->identity);
        CHECK_FALSE(first->identity == separate->identity);
        CHECK(first->byte_count == 3);
        const auto manifest = read_manifest_file(*hard_link, {3, 4, 2, 3}, plan, {}, nullptr);
        REQUIRE(manifest.has_value());
        if (const auto view = manifest->view()) {
            CHECK(view->file_info == *linked);
            CHECK(view->digest == abc_digest);
        } else {
            FAIL("Manifest has no owned view");
        }
        return read_reference_file(*primary, abc_digest, {3, 4, 2, 3}, plan, {}, nullptr);
    }();
    REQUIRE(snapshot.has_value());
    CHECK(descriptors_for(renamed).empty());
    CHECK(descriptors_for(renamed / "sub" / "value").empty());
    REQUIRE(fs::remove(renamed / "sub" / "value"));
    REQUIRE(fs::remove(renamed / "alias"));
    if (const auto view = snapshot->view()) {
        constexpr std::string_view expected = "abc";
        CHECK(std::ranges::equal(view->bytes, std::as_bytes(std::span{expected})));
        CHECK(view->digest == abc_digest);
    } else {
        FAIL("Snapshot view must survive the handles and original file");
    }
}

TEST_CASE("Linux selected directory survives replacement of its original pathname",
          "[IO4-06][FD-08][verification][property]") {
    const Tree tree;
    const auto leaf = tree.directory / "outer" / "checkpoint-1";
    const auto moved = tree.directory / "held";
    REQUIRE(fs::create_directories(leaf));
    write_file(leaf / "manifest.json", "abc");
    const auto root = LinuxReferenceRoot::open(tree.directory, {});
    REQUIRE(root);
    core::ResourcePlan plan{4096};
    const auto path = path_for("outer/checkpoint-1", plan);
    REQUIRE(path);
    const auto selected = root->open_directory(*path, plan, {});
    REQUIRE(selected);
    check_one_cloexec(leaf);

    fs::rename(leaf, moved);
    REQUIRE(fs::create_directory(leaf));
    write_file(leaf / "manifest.json", "decoy");
    const auto manifest = path_for("manifest.json", plan);
    REQUIRE(manifest);
    auto opened = selected->open_file(*manifest, plan, {});
    REQUIRE(opened);
    const auto snapshot = read_manifest_file(*opened, {3, 4, 2, 3}, plan, {}, nullptr);
    REQUIRE(snapshot);
    const auto view = snapshot->view();
    REQUIRE(view);
    CHECK(view->digest == abc_digest);
}

TEST_CASE("Linux reference lookup rejects symlinks and closes rejected native descriptors",
          "[IO4-01][verification][negative]") {
    const Tree tree;
    const auto root_path = tree.directory / "root";
    REQUIRE(fs::create_directories(root_path / "sub"));
    write_file(root_path / "sub" / "value", "abc");
    fs::create_symlink(root_path, tree.directory / "alias");
    fs::create_symlink("sub", root_path / "via");
    fs::create_symlink("sub/value", root_path / "leaf");
    REQUIRE(::mkfifo((root_path / "fifo").c_str(), S_IRUSR | S_IWUSR) == 0);
    const auto invalid_root = LinuxReferenceRoot::open(tree.directory / "alias", {});
    REQUIRE_FALSE(invalid_root.has_value());
    CHECK(invalid_root.error().code == core::FailureCode::InvalidInput);
    const auto root = LinuxReferenceRoot::open(root_path, {});
    REQUIRE(root.has_value());
    core::ResourcePlan plan{4096};
    using Reason = LinuxReferenceOpenReason;
    const std::array faults{
        std::pair<std::string_view, Reason>{"via/value", Reason::Link},
        std::pair<std::string_view, Reason>{"leaf", Reason::Link},
        std::pair<std::string_view, Reason>{"fifo", Reason::WrongType},
        std::pair<std::string_view, Reason>{"sub", Reason::WrongType},
        std::pair<std::string_view, Reason>{"missing", Reason::Open},
    };
    for (const auto& [name, expected] : faults) {
        CAPTURE(name);
        const auto path = path_for(name, plan);
        REQUIRE(path.has_value());
        std::optional<Reason> diagnostic;
        const auto rejected = root->open_file(*path, plan, {}, &diagnostic);
        REQUIRE_FALSE(rejected.has_value());
        CHECK(rejected.error().code == core::FailureCode::InvalidInput);
        CHECK(rejected.error().operation == core::Operation::PrepareModel);
        CHECK(diagnostic == expected);
    }
    CHECK(descriptors_for(root_path / "fifo").empty());
    CHECK(descriptors_for(root_path / "sub").empty());
}

TEST_CASE("Linux reference lookup admits path bytes and reports a changed handle timestamp",
          "[IO4-01][N4-02][verification][property][negative]") {
    const Tree tree;
    const auto filename = tree.directory / "name";
    write_file(filename, "abc");
    const auto root = LinuxReferenceRoot::open(tree.directory, {});
    REQUIRE(root.has_value());
    core::ResourcePlan path_plan{32};
    const auto path = path_for("name", path_plan);
    REQUIRE(path.has_value());
    const std::stop_source stop;
    REQUIRE(stop.request_stop());
    std::optional diagnostic{LinuxReferenceOpenReason::Open};
    core::ResourcePlan short_plan{4};
    const auto cancelled = root->open_file(*path, short_plan, stop.get_token(), &diagnostic);
    REQUIRE_FALSE(cancelled.has_value());
    CHECK(cancelled.error().code == core::FailureCode::Cancelled);
    CHECK_FALSE(diagnostic);
    const auto refused = root->open_file(*path, short_plan, {});
    REQUIRE_FALSE(refused.has_value());
    CHECK(refused.error().code == core::FailureCode::ResourceExceeded);
    CHECK(short_plan.allocated_bytes() == 0);
    CHECK(descriptors_for(filename).empty());

    core::ResourcePlan exact_plan{5}; // Four name bytes plus the native terminator.
    auto reader = root->open_file(*path, exact_plan, {});
    REQUIRE(reader.has_value());
    CHECK(exact_plan.allocated_bytes() == 5);
    const auto before = reader->info();
    REQUIRE(before.has_value());
    fs::last_write_time(filename, fs::last_write_time(filename) + std::chrono::seconds{17});
    const auto after = reader->info();
    REQUIRE(after.has_value());
    CHECK(after->identity == before->identity);
    CHECK(after->byte_count == before->byte_count);
    CHECK_FALSE(after->change_indicators == before->change_indicators);
}

#endif
