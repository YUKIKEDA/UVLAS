#include "case_schema_test_adapter.hpp"
#include "resolved_case_document.hpp"
#include "source_digest.hpp"

#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

#include <array>
#include <bit>
#include <charconv>
#include <cmath>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
using namespace uvlas;
using namespace uvlas::io;
using namespace uvlas::io::detail;
using Json = nlohmann::json;
constexpr std::uint64_t build_budget = std::uint64_t{1} << 50U;
constexpr std::uint64_t output_budget = 2000000;
constexpr ResolvedCaseDocumentLimits limits{1000000, output_budget, 64000000};

template <class T>
T take(core::Result<T> value) {
    if (!value) {
        throw std::runtime_error("Expected a value, code=" +
                                 std::to_string(static_cast<unsigned>(value.error().code)));
    }
    return std::move(*value);
}

struct Input {
    SiCase value;
    std::uint64_t quantities;
};

Input prepare(std::string_view wire) {
    const InputLimits syntax{
        10000000, 256,      1000000,  1000000,
        10000000, 10000000, 10000000, {ParserAllocationMode::EstimateAdmission, build_budget}};
    core::ResourcePlan plan{build_budget};
    std::optional<InputDiagnostic> diagnostic;
    const auto document = take(read_input_document(std::as_bytes(std::span{wire}), SourceId{716},
                                                   syntax, plan, diagnostic));
    std::array<std::uint64_t, static_cast<std::size_t>(SchemaQuota::Count)> schema_caps{};
    schema_caps.fill(1000000000);
    const auto checked =
        take(validate_case_schema(document, SchemaLimits{schema_caps}, plan, nullptr));
    std::array<std::uint64_t, static_cast<std::size_t>(SiQuota::Count)> si_caps{};
    si_caps.fill(1000000);
    for (const auto quota :
         {SiQuota::OwnedBytes, SiQuota::ScratchBytes, SiQuota::CumulativeBytes, SiQuota::Work}) {
        si_caps[static_cast<std::size_t>(quota)] = build_budget;
    }
    si_caps[static_cast<std::size_t>(SiQuota::UnitPower)] = 64;
    const SiLimits si_limits{si_caps,
                             {100000, 100000, 100000, 64, 64, 1000000, 2000000,
                              structure::ExactEstimateAdmission{build_budget},
                              structure::ExactTimingMode::Disabled}};
    const auto catalog = take(make_si_catalog(checked, si_limits, plan, nullptr));
    const auto count = take(catalog.size());
    std::vector<SiRequirementRow> rows;
    for (std::uint64_t i = 0; i < count; ++i) {
        const auto entry = take(catalog.entry(i));
        rows.push_back({entry.key,
                        {SiDecimalInput{"1"}, entry.recipe.dimensions},
                        SiDecimalInput{"1e-25"},
                        SiDecimalInput{"1e-25"}});
    }
    const auto request = take(make_si_request(catalog, rows, si_limits, plan, nullptr));
    return {take(convert_case_to_si(checked, catalog, request, si_limits, plan, nullptr)), count};
}

std::string sha256(std::string_view text) {
    const auto digest = source_digest(std::as_bytes(std::span{text}));
    constexpr std::string_view digits = "0123456789abcdef";
    std::string result;
    for (const auto value : digest.bytes) {
        const auto index = std::to_integer<unsigned>(value);
        result.push_back(digits[index / 16]);
        result.push_back(digits[index % 16]);
    }
    return result;
}

std::uint64_t word(const Json& value) {
    const auto& text = value.get_ref<const std::string&>();
    REQUIRE(text.size() == 16);
    std::uint64_t result = 0;
    const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), result, 16);
    REQUIRE(error == std::errc{});
    REQUIRE(end == text.data() + text.size());
    return result;
}

TEST_CASE("Resolved UTF-8 and result provenance own the same SI case", "[IO4-01][IO4-05][QA-06]") {
    using namespace uvlas::test::schema_adapter;
    auto source = root();
    source.field("units").set("length", string("mm"));
    source.field("model").field("geometry").elements[0].set("reference_length", parse("37.00"));
    const std::string wire = "\n\t" + source.wire() + "\r\n";
    const auto original = Json::parse(wire);
    constexpr std::string_view length_path = "/model/geometry/0/reference_length";
    const auto status = take(ResultStatus::completed());
    const auto bundle = [&] {
        const auto input = prepare(wire);
        const auto checked = take(input.value.source_case());
        const auto typed_root = take(checked.get().root());
        const auto geometry =
            std::get<schema::GeometryPg>(take(typed_root.model.record().geometry.at(0)).record());
        const auto typed = take(input.value.quantity(geometry.reference_length));
        const auto expected = std::get<ScaledTwoWord>(typed.enclosure);
        core::ResourcePlan plan{output_budget};
        auto documents = take(
            write_resolved_result_documents(input.value, status, limits, output_budget, plan, {}));
        const auto bytes = documents.resolved_case.text();
        const auto resolved = Json::parse(bytes);
        const auto& rows = resolved.at("quantities");
        const auto& defaults = resolved.at("structural_defaults");
        CHECK(resolved.at("original_units").at("length") == "mm");
        CHECK(defaults.at("provenance_missing") == true);
        CHECK(defaults.at("initial_missing") == true);
        CHECK(defaults.at("restart_missing") == true);
        const Json empty_provenance = {
            {"description", ""}, {"references", Json::array()}, {"assumptions", Json::array()}};
        CHECK(defaults.at("applied_provenance") == empty_provenance);
        REQUIRE(rows.size() == input.quantities);
        std::set<std::string> unique_paths;
        const Json* length = nullptr;
        for (std::size_t i = 0; i < rows.size(); ++i) {
            const auto& row = rows.at(i);
            const auto& path = row.at("path").get_ref<const std::string&>();
            CHECK(unique_paths.insert(path).second);
            CHECK(original.contains(Json::json_pointer{path}));
            CHECK(row.at("ledger_origin") == std::to_string(i));
            if (path == length_path) {
                length = &row;
            }
        }
        REQUIRE(length);
        if (length) {
            CHECK(length->at("dimension") == Json::array({1, 0, 0, 0}));
            CHECK(length->at("original").at("decimal").at("digits") == "37");
            CHECK(length->at("original").at("token") == "37.00");
            const Json requests = {{"si_scale", "1"}, {"absolute", "1e-25"}, {"relative", "1e-25"}};
            CHECK(length->at("requirements") == requests);
            const auto& enclosure = length->at("enclosure");
            CHECK(enclosure.at("representation") == "TwoWord");
            CHECK(enclosure.at("exponent2") == std::to_string(expected.exponent2));
            CHECK(word(enclosure.at("high_bits")) == std::bit_cast<std::uint64_t>(expected.high));
            CHECK(word(enclosure.at("low_bits")) == std::bit_cast<std::uint64_t>(expected.low));
            CHECK(word(enclosure.at("radius_bits")) ==
                  std::bit_cast<std::uint64_t>(expected.radius));
            CHECK(expected.low != 0.0);
            CHECK(std::ldexp(expected.high, static_cast<int>(expected.exponent2)) == 0.037);
        }
        CHECK(plan.allocated_bytes() == bytes.size() + documents.result.status().size() +
                                            documents.result.provenance().size());
        const ResolvedCaseDocumentLimits tight{input.quantities, bytes.size(), limits.maximum_work};
        core::ResourcePlan exact{bytes.size()};
        const auto replay = take(write_resolved_case_document(input.value, tight, exact, {}));
        CHECK(replay.text() == bytes);
        CHECK(exact.allocated_bytes() == bytes.size());

        const auto reject = [&](ResolvedCaseDocumentLimits caps, std::uint64_t memory,
                                core::FailureCode code, const std::stop_token& stop = {}) {
            core::ResourcePlan denied{memory};
            const auto rejected = write_resolved_case_document(input.value, caps, denied, stop);
            REQUIRE_FALSE(rejected);
            CHECK(rejected.error().code == code);
            CHECK(rejected.error().operation == core::Operation::WriteResults);
            CHECK(denied.allocated_bytes() == 0);
        };
        const auto resource = core::FailureCode::ResourceExceeded;
        reject({input.quantities - 1, output_budget, limits.maximum_work}, output_budget, resource);
        reject({input.quantities, bytes.size() - 1, limits.maximum_work}, output_budget, resource);
        reject({input.quantities, bytes.size(), bytes.size()}, output_budget, resource);
        reject(tight, bytes.size() - 1, resource);
        // NOLINTNEXTLINE(misc-const-correctness): MSVC request_stop() is a non-const member.
        std::stop_source cancelled;
        REQUIRE(cancelled.request_stop());
        reject(tight, output_budget, core::FailureCode::Cancelled, cancelled.get_token());
        core::ResourcePlan late_failure{output_budget};
        const auto rejected =
            write_resolved_result_documents(input.value, status, tight, 0, late_failure, {});
        REQUIRE_FALSE(rejected);
        CHECK(rejected.error().code == resource);
        CHECK(late_failure.allocated_bytes() == bytes.size());
        return documents;
    }();
    const auto resolved = Json::parse(bundle.resolved_case.text());
    const auto result = Json::parse(bundle.result.status());
    const auto provenance = Json::parse(bundle.result.provenance());
    CHECK(resolved.at("original_case_utf8") == wire);
    CHECK(resolved.at("source").at("sha256") == sha256(wire));
    CHECK(result.at("completed") == true);
    CHECK(result.at("provenance_sha256") == sha256(bundle.result.provenance()));
    CHECK(provenance.at("original_case").at("source_id") == "716");
    CHECK(provenance.at("original_case").at("byte_count") == std::to_string(wire.size()));
    CHECK(provenance.at("original_case").at("sha256") == sha256(wire));
    CHECK(provenance.at("resolved_case_sha256") == sha256(bundle.resolved_case.text()));
}
} // namespace
