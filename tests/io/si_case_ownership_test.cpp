#include "case_schema_test_adapter.hpp"
#include "si_case_ownership_snapshot.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <uvlas/io/si_case.hpp>
#include <vector>
namespace {
using namespace uvlas;
using namespace uvlas::test::schema_adapter;
constexpr std::uint64_t budget = std::uint64_t{1} << 50U;
template <class T>
T take(core::Result<T> value) {
    if (!value)
        throw std::runtime_error("SI unit control rejected, code=" +
                                 std::to_string(static_cast<unsigned>(value.error().code)));
    return std::move(*value);
}
template <class T>
T si_phase(core::Result<T> value, std::string_view phase, const io::SiDiagnostic& diag) {
    if (!value) {
        std::string message{phase};
        message += " code=" + std::to_string(static_cast<unsigned>(value.error().code));
        const auto detail = diag.detail();
        if (detail)
            message += " reason=" + std::to_string(static_cast<unsigned>(detail->reason));
        const auto path = diag.path();
        if (path)
            message += " path=" + std::string(*path);
        throw std::runtime_error(message);
    }
    return std::move(*value);
}
io::SiLimits unit_limits() {
    std::array<std::uint64_t, static_cast<std::size_t>(io::SiQuota::Count)> q{};
    q.fill(1000000);
    for (const auto quota : {io::SiQuota::OwnedBytes, io::SiQuota::ScratchBytes,
                             io::SiQuota::CumulativeBytes, io::SiQuota::Work})
        q[static_cast<std::size_t>(quota)] = budget;
    q[static_cast<std::size_t>(io::SiQuota::UnitPower)] = 64;
    return {q,
            {100000, 100000, 100000, 64, 64, 1000000, 2000000,
             structure::ExactEstimateAdmission{budget}, structure::ExactTimingMode::Disabled}};
}
io::SchemaCase checked(const JsonToken& source) {
    const auto wire = source.wire();
    const io::InputLimits syntax{
        10000000, 256,      1000000,  1000000,
        10000000, 10000000, 10000000, {io::ParserAllocationMode::EstimateAdmission, budget}};
    core::ResourcePlan plan{budget};
    std::optional<io::InputDiagnostic> diag;
    const auto doc = take(io::read_input_document(
        std::as_bytes(std::span{wire.data(), wire.size()}), io::SourceId{716}, syntax, plan, diag));
    std::array<std::uint64_t, static_cast<std::size_t>(io::SchemaQuota::Count)> caps{};
    caps.fill(1000000000);
    io::SchemaDiagnostic schema;
    auto result = io::validate_case_schema(doc, io::SchemaLimits{caps}, plan, &schema);
    if (!result) {
        std::string message =
            "schema code=" + std::to_string(static_cast<unsigned>(result.error().code));
        const auto detail = schema.detail();
        if (detail)
            message += " reason=" + std::to_string(static_cast<unsigned>(detail->reason));
        const auto path = schema.path();
        if (path)
            message += " path=" + std::string(*path);
        throw std::runtime_error(message);
    }
    return std::move(*result);
}
std::vector<io::SiRequirementRow> requirements(const io::SiCatalog& catalog,
                                               std::optional<io::SiQuantityKey> exact = {}) {
    std::vector<io::SiRequirementRow> rows;
    for (std::uint64_t i = 0; i < take(catalog.size()); ++i) {
        const auto entry = take(catalog.entry(i));
        const bool zero = exact && *exact == entry.key;
        rows.push_back({entry.key,
                        {io::SiDecimalInput{"1"}, entry.recipe.dimensions},
                        io::SiDecimalInput{zero ? "0" : "1e-25"},
                        io::SiDecimalInput{zero ? "0" : "1e-25"}});
    }
    return rows;
}

JsonToken ownership_source() {
    auto source = root();
    source.field("units").set("length", string("mm"));
    source.field("model").field("geometry").elements[0].set("reference_length", parse("37"));
    return source;
}
struct Owned {
    io::SiCase value;
    std::uint64_t count;
};
Owned prepare_owned() {
    const auto schema = checked(ownership_source());
    core::ResourcePlan plan{budget};
    const auto limits = unit_limits();
    io::SiDiagnostic diagnostic;
    const auto catalog =
        si_phase(io::make_si_catalog(schema, limits, plan, &diagnostic), "catalog", diagnostic);
    const auto count = take(catalog.size());
    const auto rows = requirements(catalog);
    const auto request = si_phase(io::make_si_request(catalog, rows, limits, plan, &diagnostic),
                                  "request", diagnostic);
    const auto before = plan.allocated_bytes();
    auto value =
        si_phase(io::convert_case_to_si(schema, catalog, request, limits, plan, &diagnostic),
                 "convert", diagnostic);
    const auto report = take(value.memory_report());
    CHECK(plan.allocated_bytes() - before == report.cumulative_bytes);
    CHECK(report.owned_bytes <= report.cumulative_bytes);
    return {std::move(value), count};
}
auto snapshot(const Owned& value) {
    const auto result = uvlas::test::si_ownership_snapshot::serialize(value.value);
    REQUIRE(result.keys.size() == value.count);
    return result;
}
} // namespace
TEST_CASE("SI results retain every public payload after all producing owners and plan die",
          "[IO4-01][property][verification]") {
    auto original = prepare_owned();
    const auto before = snapshot(original);
    const auto source = take(original.value.source_case());
    CHECK(take(source.get().root()).units.record().length.text() == "mm");
    auto moved = std::move(original.value);
    const auto absent = original.value.source_case();
    REQUIRE_FALSE(absent);
    CHECK(absent.error().code == core::FailureCode::InvalidInput);
    const auto moved_snapshot = uvlas::test::si_ownership_snapshot::serialize(moved);
    CHECK(moved_snapshot == before);
    CHECK(moved_snapshot.keys.size() == original.count);
    original.value = std::move(moved);
    // NOLINTNEXTLINE(bugprone-use-after-move): the public moved-from query must reject safely.
    REQUIRE_FALSE(moved.source_case());
    CHECK(snapshot(original) == before);
}
TEST_CASE("SI quantity access rejects a foreign typed owner even with equal source bytes",
          "[IO4-01][negative][property]") {
    const auto original = prepare_owned();
    const auto before = snapshot(original);
    const auto foreign = checked(ownership_source());
    const auto own_source = take(original.value.source_case());
    CHECK(take(own_source.get().source_identity()).digest ==
          take(foreign.source_identity()).digest);
    const auto own_root = take(own_source.get().root());
    const auto foreign_root = take(foreign.root());
    const auto own_geometry =
        std::get<io::schema::GeometryPg>(take(own_root.model.record().geometry.at(0)).record());
    const auto foreign_geometry =
        std::get<io::schema::GeometryPg>(take(foreign_root.model.record().geometry.at(0)).record());
    CHECK(original.value.quantity(own_geometry.reference_length));
    const auto rejected = original.value.quantity(foreign_geometry.reference_length);
    REQUIRE_FALSE(rejected);
    CHECK(rejected.error().code == core::FailureCode::InvalidInput);
    CHECK(snapshot(original) == before);
}
TEST_CASE("SI requests reject keys from another catalog of the identical checked case",
          "[IO4-01][negative][property]") {
    const auto original = prepare_owned();
    const auto before = snapshot(original);
    const auto source = checked(ownership_source());
    const auto limits = unit_limits();
    core::ResourcePlan plan{budget};
    io::SiDiagnostic diag;
    // Retain the foreign keys while releasing their catalog, then issue a fresh
    // catalog. A retained owned generation must prevent address-reuse/ABA acceptance.
    const auto foreign_rows = [&] {
        const auto catalog =
            si_phase(io::make_si_catalog(source, limits, plan, &diag), "foreign catalog", diag);
        return requirements(catalog);
    }();
    const auto catalog =
        si_phase(io::make_si_catalog(source, limits, plan, &diag), "catalog", diag);
    const auto good = requirements(catalog);
    REQUIRE(foreign_rows.size() == good.size());
    const auto rejected = io::make_si_request(catalog, foreign_rows, limits, plan, &diag);
    REQUIRE_FALSE(rejected);
    CHECK(rejected.error().code == core::FailureCode::InvalidInput);
    const auto accepted = io::make_si_request(catalog, good, limits, plan, &diag);
    REQUIRE(accepted);
    CHECK(snapshot(original) == before);
}
TEST_CASE("SI failed calls preserve every old value and successful allocations are not refunded",
          "[IO4-01][negative][property]") {
    const auto original = prepare_owned();
    const auto before = snapshot(original);
    const auto source = checked(ownership_source());
    const auto source_identity = take(source.source_identity());
    const auto limits = unit_limits();
    core::ResourcePlan plan{budget};
    io::SiDiagnostic diag;
    const auto catalog =
        si_phase(io::make_si_catalog(source, limits, plan, &diag), "catalog", diag);
    const auto correct = requirements(catalog);
    REQUIRE(correct.size() > 1);
    for (const auto defect : {0, 1, 2}) {
        INFO(defect);
        auto rows = correct;
        if (defect == 0)
            rows.pop_back();
        if (defect == 1)
            rows.back() = rows.front();
        if (defect == 2)
            rows.front().scale.dimensions.exponents[0] = 17;
        const auto charged = plan.allocated_bytes();
        const auto rejected = io::make_si_request(catalog, rows, limits, plan, &diag);
        REQUIRE_FALSE(rejected);
        CHECK(rejected.error().code == core::FailureCode::InvalidInput);
        CHECK(plan.allocated_bytes() >= charged);
        CHECK(snapshot(original) == before);
        CHECK(take(source.source_identity()).digest == source_identity.digest);
    }
    const auto before_request = plan.allocated_bytes();
    std::optional<io::SiRequest> request;
    request.emplace(
        si_phase(io::make_si_request(catalog, correct, limits, plan, &diag), "request", diag));
    const auto after_request = plan.allocated_bytes();
    CHECK(after_request - before_request == take(request->memory_report()).cumulative_bytes);
    core::ResourcePlan empty{0};
    const auto failed = io::convert_case_to_si(source, catalog, *request, limits, empty, &diag);
    REQUIRE_FALSE(failed);
    CHECK(failed.error().code == core::FailureCode::ResourceExceeded);
    CHECK(empty.allocated_bytes() == 0);
    CHECK(snapshot(original) == before);
    request.reset();
    CHECK(plan.allocated_bytes() == after_request);
    CHECK(snapshot(original) == before);
}

TEST_CASE("SI host request lexical boundaries normalize before exponent admission",
          "[IO4-01][verification][negative][property]") {
    const auto source = checked(ownership_source());
    const auto base = unit_limits();
    std::array<std::uint64_t, static_cast<std::size_t>(io::SiQuota::Count)> quotas{};
    for (std::size_t i = 0; i < quotas.size(); ++i)
        quotas[i] = take(base.get(static_cast<io::SiQuota>(i)));
    quotas[static_cast<std::size_t>(io::SiQuota::ExponentMagnitude)] = 64;
    const io::SiLimits limits{quotas, base.exact()};
    core::ResourcePlan plan{budget};
    io::SiDiagnostic diagnostic;
    const auto catalog =
        si_phase(io::make_si_catalog(source, limits, plan, &diagnostic), "catalog", diagnostic);
    const auto good = requirements(catalog);
    REQUIRE_FALSE(good.empty());
    const auto normalized_one = std::string{"1"} + std::string(100, '0') + "e-100";
    {
        auto rows = good;
        rows.front().scale.value = io::SiDecimalInput{normalized_one};
        rows.front().absolute = io::SiDecimalInput{"1"};
        rows.front().relative = io::SiDecimalInput{"1"};
        const auto accepted = io::make_si_request(catalog, rows, limits, plan, &diagnostic);
        REQUIRE(accepted);
    }
    {
        auto rows = good;
        rows.front().absolute = io::SiDecimalInput{"-0e999999999999999999999999"};
        rows.front().relative = io::SiDecimalInput{"-0e999999999999999999999999"};
        const auto accepted = io::make_si_request(catalog, rows, limits, plan, &diagnostic);
        REQUIRE(accepted);
    }
    {
        auto rows = good;
        rows.front().relative = io::SiDecimalInput{"-1e-9223372036854775807"};
        const auto rejected = io::make_si_request(catalog, rows, limits, plan, &diagnostic);
        REQUIRE_FALSE(rejected);
        CHECK(rejected.error().code == core::FailureCode::InvalidInput);
    }
    for (const auto bad : {io::SiQuota::Count, static_cast<io::SiQuota>(255)}) {
        const auto rejected = limits.get(bad);
        REQUIRE_FALSE(rejected);
        CHECK(rejected.error().code == core::FailureCode::InvalidInput);
    }
}

TEST_CASE("SI unproven external memory policy rejects all three factories without reservation",
          "[IO4-01][verification][negative][property]") {
    const auto source = checked(ownership_source());
    const auto normal = unit_limits();
    core::ResourcePlan plan{budget};
    io::SiDiagnostic diagnostic;
    const auto catalog =
        si_phase(io::make_si_catalog(source, normal, plan, &diagnostic), "catalog", diagnostic);
    const auto rows = requirements(catalog);
    const auto request = si_phase(io::make_si_request(catalog, rows, normal, plan, &diagnostic),
                                  "request", diagnostic);
    std::array<std::uint64_t, static_cast<std::size_t>(io::SiQuota::Count)> quotas{};
    for (std::size_t i = 0; i < quotas.size(); ++i)
        quotas[i] = take(normal.get(static_cast<io::SiQuota>(i)));
    auto exact = normal.exact();
    exact.external_memory = structure::ExactRequireProvenBound{};
    const io::SiLimits proven{quotas, exact};
    const auto reject = [&](auto call) {
        const auto charged = plan.allocated_bytes(), remaining = plan.remaining_bytes();
        const auto result = call();
        REQUIRE_FALSE(result);
        CHECK(result.error().code == core::FailureCode::Unsupported);
        CHECK(plan.allocated_bytes() == charged);
        CHECK(plan.remaining_bytes() == remaining);
    };
    reject([&] { return io::make_si_catalog(source, proven, plan, &diagnostic); });
    reject([&] { return io::make_si_request(catalog, rows, proven, plan, &diagnostic); });
    reject([&] {
        return io::convert_case_to_si(source, catalog, request, proven, plan, &diagnostic);
    });
}
