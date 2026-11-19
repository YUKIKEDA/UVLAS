#include "../support/allocation_audit.hpp"
#include "reference_path.hpp"

#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <optional>
#include <ranges>
#include <stop_token>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {
using namespace uvlas;
using namespace uvlas::io::detail;
using namespace uvlas::test_support;
using Reason = ReferencePathReason;

constexpr ReferencePathLimits ordinary_limits{1024, 1024, 64};

[[nodiscard]] std::optional<std::string> stack_reference(std::string_view input) {
    std::vector<std::string_view> retained;
    for (const auto part : input | std::views::split('/')) {
        const std::string_view component(part.begin(), part.end());
        if (component == "..") {
            if (retained.empty()) {
                return std::nullopt;
            }
            retained.pop_back();
        } else if (!component.empty() && component != ".") {
            retained.push_back(component);
        }
    }
    if (retained.empty()) {
        return std::nullopt;
    }
    std::string output;
    for (const auto component : retained) {
        if (!output.empty()) {
            output += '/';
        }
        output.append(component);
    }
    return output;
}

} // namespace

TEST_CASE("Reference path normalization matches the frozen examples",
          "[CS-04][IO4-01][verification][negative]") {
    std::ifstream fixture(UVLAS_REFERENCE_PATH_FIXTURE, std::ios::binary);
    REQUIRE(fixture.is_open());
    const auto examples = nlohmann::json::parse(fixture).at("paths");
    REQUIRE(examples.size() == 11);
    constexpr std::array error_names{
        std::pair{"EscapeRoot", Reason::EscapeRoot},
        std::pair{"RootedPath", Reason::RootedPath},
        std::pair{"StreamPath", Reason::StreamPath},
        std::pair{"MissingFile", Reason::MissingFile},
        std::pair{"InvalidPath", Reason::InvalidPath},
    };
    const auto stop = std::stop_token{};
    for (const auto& example : examples) {
        const auto path = example.at("declared").get<std::string>();
        CAPTURE(path);
        core::ResourcePlan plan{1024};
        std::optional<ReferencePathDiagnostic> diagnostic;
        const auto result =
            normalize_reference_path(path, ordinary_limits, plan, stop, &diagnostic);
        if (example.contains("normalized")) {
            const auto expected = example.at("normalized").get<std::string>();
            REQUIRE(result.has_value());
            CHECK(result->text() == expected);
            CHECK(plan.allocated_bytes() == expected.size());
            CHECK_FALSE(diagnostic.has_value());
        } else {
            REQUIRE_FALSE(result.has_value());
            CHECK(result.error().code == core::FailureCode::InvalidInput);
            CHECK(result.error().operation == core::Operation::PrepareModel);
            const auto expected = example.at("failure").get<std::string>();
            const auto entry = std::ranges::find_if(
                error_names, [&](auto pair) { return pair.first == expected; });
            REQUIRE(entry != error_names.end());
            if (diagnostic) {
                CHECK(diagnostic->reason == entry->second);
            } else {
                FAIL("Missing lexical diagnostic");
            }
            CHECK(plan.allocated_bytes() == 0);
        }
    }
}

TEST_CASE("Reference path reverse scan agrees with an independent stack",
          "[CS-04][IO4-01][property]") {
    constexpr std::array<std::string_view, 5> components{"alpha", "B ", "..", ".", "z"};
    const auto stop = std::stop_token{};
    for (const auto first : components) {
        for (const auto second : components) {
            for (const auto third : components) {
                for (const auto fourth : components) {
                    std::string path(first);
                    path.append("/").append(second).append("//").append(third).append("/").append(
                        fourth);
                    const auto expected = stack_reference(path);
                    core::ResourcePlan plan{1024};
                    const auto result =
                        normalize_reference_path(path, ordinary_limits, plan, stop, nullptr);
                    CAPTURE(path);
                    REQUIRE(result.has_value() == expected.has_value());
                    if (expected) {
                        CHECK(result->text() == *expected);
                    }
                }
            }
        }
    }
}

TEST_CASE("Reference paths reserve only the exact owned payload",
          "[CS-04][N4-02][N4-09][negative][property]") {
    const auto stop = std::stop_token{};
    constexpr std::string_view input = "drop/../A//B ";
    struct Budget {
        ReferencePathLimits limits;
        std::uint64_t plan_size;
        Reason reason;
    };
    constexpr std::array shortages{
        Budget{{input.size() - 1, 4, 2}, 4, Reason::InputBytes},
        Budget{{input.size(), 3, 2}, 4, Reason::NormalizedBytes},
        Budget{{input.size(), 4, 1}, 4, Reason::NormalizedDepth},
        Budget{{input.size(), 4, 2}, 3, Reason::PlanBytes},
    };
    for (const auto& shortage : shortages) {
        core::ResourcePlan plan{shortage.plan_size};
        std::optional<ReferencePathDiagnostic> diagnostic;
        AllocationAudit audit{cpp_backend, AuditThreadDomain::AllThreads};
        const auto result =
            normalize_reference_path(input, shortage.limits, plan, stop, &diagnostic);
        const auto allocation = audit.finish();
        CHECK(allocation.status == AllocationAuditStatus::Clean);
        CHECK(allocation.cpp_allocation_attempts == 0);
        REQUIRE_FALSE(result.has_value());
        CHECK(result.error().code == core::FailureCode::ResourceExceeded);
        if (diagnostic) {
            CHECK(diagnostic->reason == shortage.reason);
        } else {
            FAIL("Missing resource diagnostic");
        }
        CHECK(plan.allocated_bytes() == 0);
    }

    core::ResourcePlan plan{4};
    AllocationAudit audit{cpp_backend, AuditThreadDomain::AllThreads};
    const auto result = normalize_reference_path(input, {input.size(), 4, 2}, plan, stop, nullptr);
    const auto allocation = audit.finish();
    CHECK(allocation.status == AllocationAuditStatus::AllocationObserved);
    CHECK(allocation.cpp_allocation_attempts == 1);
    REQUIRE(result.has_value());
    CHECK(result->text() == "A/B ");
    CHECK(plan.allocated_bytes() == 4);
}

TEST_CASE("Reference path failure points and lifetime are stable",
          "[CS-04][IO4-01][negative][property]") {
    const auto stop = std::stop_token{};
    const auto obtain = [&] {
        const std::string input = " Mixed /Decomposed-e\xCC\x81";
        core::ResourcePlan plan{1024};
        return normalize_reference_path(input, ordinary_limits, plan, stop, nullptr);
    };
    auto saved = obtain();
    REQUIRE(saved.has_value());
    auto moved = std::move(*saved);
    CHECK(moved.text() == " Mixed /Decomposed-e\xCC\x81");

    core::ResourcePlan plan{1024};
    std::optional<ReferencePathDiagnostic> diagnostic;
    const auto escaped =
        normalize_reference_path("a/../b/../../c", ordinary_limits, plan, stop, &diagnostic);
    REQUIRE_FALSE(escaped.has_value());
    if (diagnostic) {
        CHECK(diagnostic->reason == Reason::EscapeRoot);
        CHECK(diagnostic->input_offset == 10);
    } else {
        FAIL("Missing escape diagnostic");
    }

    // NOLINTNEXTLINE(misc-const-correctness): MSVC 14.51 request_stop() is non-const.
    std::stop_source source;
    REQUIRE(source.request_stop());
    const auto cancelled =
        normalize_reference_path("valid", ordinary_limits, plan, source.get_token(), &diagnostic);
    REQUIRE_FALSE(cancelled.has_value());
    CHECK(cancelled.error().code == core::FailureCode::Cancelled);
    CHECK_FALSE(diagnostic.has_value());
    CHECK(plan.allocated_bytes() == 0);
    CHECK(moved.text() == " Mixed /Decomposed-e\xCC\x81");
}
