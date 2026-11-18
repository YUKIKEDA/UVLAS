#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cstdint>
#include <limits>
#include <optional>
#include <string_view>
#include <uvlas/io/case_schema.hpp>

namespace {
uvlas::io::SchemaLimits schema_limits() {
    std::array<std::uint64_t, static_cast<std::size_t>(uvlas::io::SchemaQuota::Count)> values{};
    values.fill(1000000);
    values[static_cast<std::size_t>(uvlas::io::SchemaQuota::TypedBytes)] = 64ULL * 1024ULL * 1024;
    values[static_cast<std::size_t>(uvlas::io::SchemaQuota::ScratchBytes)] = 16ULL * 1024ULL * 1024;
    values[static_cast<std::size_t>(uvlas::io::SchemaQuota::Work)] = 1000000000;
    values[static_cast<std::size_t>(uvlas::io::SchemaQuota::ExactDigits)] = 4096;
    return uvlas::io::SchemaLimits{values};
}
} // namespace
TEST_CASE("Schema quotas reject sentinel and unknown enumerators", "[unit][negative][IO4-01]") {
    const auto limits = schema_limits();
    CHECK(limits.get(uvlas::io::SchemaQuota::TypedBytes).has_value());
    const auto sentinel = limits.get(uvlas::io::SchemaQuota::Count);
    REQUIRE_FALSE(sentinel.has_value());
    CHECK(sentinel.error().code == uvlas::core::FailureCode::InvalidInput);
    const auto unknown = limits.get(static_cast<uvlas::io::SchemaQuota>(255));
    REQUIRE_FALSE(unknown.has_value());
    CHECK(unknown.error().code == uvlas::core::FailureCode::InvalidInput);
}
TEST_CASE("Nested null values cannot produce a checked SchemaCase", "[negative][IO4-01]") {
    using namespace uvlas;
    constexpr std::string_view source =
        R"({"schema":"uvlas.case/1","units":null,"model":null,"flow":null,"analysis":null,"observations":null,"numerics":null,"resources":null,"output":null})";
    const io::InputLimits input_limits{
        10000,
        64,
        1000,
        1000,
        1000,
        1000,
        10000,
        {io::ParserAllocationMode::EstimateAdmission, std::numeric_limits<std::uint64_t>::max()}};
    core::ResourcePlan input_plan{1024ULL * 1024};
    std::optional<io::InputDiagnostic> input_diagnostic;
    auto document =
        io::read_input_document(std::as_bytes(std::span{source.data(), source.size()}),
                                io::SourceId{0}, input_limits, input_plan, input_diagnostic);
    if (!document) {
        FAIL("Envelope control must succeed");
        return;
    }
    core::ResourcePlan plan{64ULL * 1024ULL * 1024};
    io::SchemaDiagnostic diagnostic;
    const auto checked = io::validate_case_schema(*document, schema_limits(), plan, &diagnostic);
    REQUIRE_FALSE(checked.has_value());
    CHECK(checked.error().code == core::FailureCode::InvalidInput);
    const auto detail = diagnostic.detail();
    if (!detail) {
        FAIL("Schema failure needs owned detail");
        return;
    }
    CHECK(detail->reason == io::SchemaReason::WrongKind);
    CHECK(document->view().has_value());
}

#include "case_schema_corpus.hpp"
#include "case_schema_payload.hpp"
#include "case_schema_test_adapter.hpp"
namespace {
uvlas::io::detail::FieldType declared_type(std::string_view spelling) {
    const auto types = uvlas::io::detail::schema_types();
    for (std::size_t i = 0; i < types.size(); ++i)
        if (types[i].spelling == spelling)
            return static_cast<uvlas::io::detail::FieldType>(i);
    throw std::invalid_argument("Unknown declared schema type");
}
uvlas::core::Result<uvlas::io::detail::SchemaPayload>
scoped_payload(std::string_view value, std::string_view type, uvlas::io::SchemaLimits limits,
               uvlas::core::ResourcePlan& plan, uvlas::io::SchemaDiagnostic* diagnostic = nullptr) {
    using namespace uvlas;
    const auto wire = std::string{"{\"schema\":\"uvlas.case/1\",\"units\":"} + std::string{value} +
                      ",\"model\":null,\"flow\":null,\"analysis\":null,\"observations\":null,"
                      "\"numerics\":null,\"resources\":null,\"output\":null}";
    const io::InputLimits syntax{
        1000000,
        256,
        100000,
        100000,
        1000000,
        1000000,
        1000000,
        {io::ParserAllocationMode::EstimateAdmission, std::numeric_limits<std::uint64_t>::max()}};
    core::ResourcePlan source_plan{64ULL * 1024ULL * 1024};
    std::optional<io::InputDiagnostic> source_diagnostic;
    auto document =
        io::read_input_document(std::as_bytes(std::span{wire.data(), wire.size()}),
                                io::SourceId{21}, syntax, source_plan, source_diagnostic);
    if (!document)
        return std::unexpected(document.error());
    const auto view = document->view();
    if (!view)
        return std::unexpected(view.error());
    const auto object = view->root.object();
    if (!object)
        return std::unexpected(
            core::Failure{core::FailureCode::InvalidInput, core::Operation::PrepareModel});
    const auto payload = object->find("units");
    if (!payload)
        return std::unexpected(
            core::Failure{core::FailureCode::InvalidInput, core::Operation::PrepareModel});
    return io::detail::validate_schema_payload(
        *document, {*payload, declared_type(type), 1, io::SchemaDimension{{0, 0, 0, 0}}, false},
        limits, plan, diagnostic);
}
uvlas::io::SchemaLimits with_quota(uvlas::io::SchemaQuota quota, std::uint64_t value) {
    std::array<std::uint64_t, static_cast<std::size_t>(uvlas::io::SchemaQuota::Count)> values{};
    const auto original = schema_limits();
    for (std::size_t i = 0; i < values.size(); ++i)
        values[i] = *original.get(static_cast<uvlas::io::SchemaQuota>(i));
    values[static_cast<std::size_t>(quota)] = value;
    return uvlas::io::SchemaLimits{values};
}
} // namespace
TEST_CASE("The shared exact kernel subtracts the smallest signed decimal exponent safely",
          "[verification][negative][IO4-01]") {
    using namespace uvlas;
    constexpr auto value =
        R"({"start":1.1e-9223372036854775807,"end":2.2e-9223372036854775807,"step":1.1e-9223372036854775807,"restart":"initial"})";
    core::ResourcePlan plan{64ULL * 1024ULL * 1024};
    io::SchemaDiagnostic diagnostic;
    auto result = scoped_payload(value, "TimeInterval",
                                 with_quota(io::SchemaQuota::ExponentMagnitude, 1ULL << 63), plan,
                                 &diagnostic);
    if (const auto detail = diagnostic.detail()) {
        INFO(static_cast<unsigned>(detail->reason));
        if (const auto path = diagnostic.path())
            INFO(*path);
    }
    REQUIRE(result.has_value());
    const auto record = io::detail::RecordAccess<io::detail::FieldType::T409>::read(
        result->data.values().data(), result->root);
    CHECK(record.start.number().exponent10 == std::numeric_limits<std::int64_t>::min());
    CHECK(record.start.number().digits == "11");
}
TEST_CASE("Positive polynomial table and expression charges have exact independent limits",
          "[verification][negative][IO4-01]") {
    using namespace uvlas;
    struct Control {
        std::string_view value, type;
        io::SchemaQuota quota;
        std::uint64_t needed;
    };
    const std::array controls{
        Control{
            R"({"kind":"polynomial","coordinates":"chart","origin":[0,0,0],"lengths":[1,1,1],"components":[[{"powers":[1,0,0],"coefficient":1}]],"error":[0]})",
            "ReferenceField", io::SchemaQuota::PolynomialTerms, 1},
        Control{R"({"kind":"pchip","times":[0,1],"values":[[0],[0]],"outside":{"kind":"reject"}})",
                "TimeSeries", io::SchemaQuota::TableSamples, 2},
        Control{
            R"({"domain":{"kind":"reference","coordinates":"chart","lower":[0,0,0],"upper":[1,1,1]},"nodes":[{"op":"constant","value":1,"dimensions":[0,0,0,0]},{"op":"integer_power","operand":0,"exponent":3}],"outputs":[1]})",
            "Expression", io::SchemaQuota::ExpressionNodes, 2},
        Control{
            R"({"domain":{"kind":"reference","coordinates":"chart","lower":[0,0,0],"upper":[1,1,1]},"nodes":[{"op":"constant","value":1,"dimensions":[0,0,0,0]},{"op":"integer_power","operand":0,"exponent":3}],"outputs":[1]})",
            "Expression", io::SchemaQuota::ExpressionEdges, 1},
        Control{
            R"({"domain":{"kind":"reference","coordinates":"chart","lower":[0,0,0],"upper":[1,1,1]},"nodes":[{"op":"constant","value":1,"dimensions":[0,0,0,0]},{"op":"integer_power","operand":0,"exponent":3}],"outputs":[1]})",
            "Expression", io::SchemaQuota::ExpressionDepth, 2},
        Control{
            R"({"domain":{"kind":"reference","coordinates":"chart","lower":[0,0,0],"upper":[1,1,1]},"nodes":[{"op":"constant","value":1,"dimensions":[0,0,0,0]},{"op":"integer_power","operand":0,"exponent":3}],"outputs":[1]})",
            "Expression", io::SchemaQuota::ExpressionPower, 3},
    };
    for (const auto& control : controls) {
        INFO(control.type);
        INFO(static_cast<unsigned>(control.quota));
        core::ResourcePlan equal_plan{64ULL * 1024ULL * 1024};
        io::SchemaDiagnostic equal_diagnostic;
        auto equal =
            scoped_payload(control.value, control.type, with_quota(control.quota, control.needed),
                           equal_plan, &equal_diagnostic);
        if (const auto detail = equal_diagnostic.detail()) {
            INFO("reason " << static_cast<unsigned>(detail->reason));
            if (const auto path = equal_diagnostic.path())
                INFO(*path);
        }
        REQUIRE(equal.has_value());
        core::ResourcePlan short_plan{64ULL * 1024ULL * 1024};
        io::SchemaDiagnostic short_diagnostic;
        const auto short_result = scoped_payload(control.value, control.type,
                                                 with_quota(control.quota, control.needed - 1),
                                                 short_plan, &short_diagnostic);
        REQUIRE_FALSE(short_result.has_value());
        CHECK(short_result.error().code == core::FailureCode::ResourceExceeded);
        const auto detail = short_diagnostic.detail();
        if (!detail) {
            FAIL_CHECK("Quota rejection needs detail");
            continue;
        }
        CHECK(detail->quota == control.quota);
        CHECK(detail->required == control.needed);
    }
}
TEST_CASE("Exact decimal intervals retain a positive step beyond binary64 resolution",
          "[verification][IO4-01]") {
    using namespace uvlas;
    constexpr auto value =
        R"({"start":1,"end":1.0000000000000000000000001,"step":0.0000000000000000000000001,"restart":"initial"})";
    core::ResourcePlan plan{64ULL * 1024 * 1024};
    io::SchemaDiagnostic diagnostic;
    auto result = scoped_payload(value, "TimeInterval", schema_limits(), plan, &diagnostic);
    REQUIRE(result.has_value());
    const auto record = io::detail::RecordAccess<io::detail::FieldType::T409>::read(
        result->data.values().data(), result->root);
    CHECK(record.start.number().digits == "1");
    CHECK(record.start.number().exponent10 == 0);
    CHECK(record.end.number().digits == "10000000000000000000000001");
    CHECK(record.end.number().exponent10 == -25);
    CHECK(record.step.number().digits == "1");
    CHECK(record.step.number().exponent10 == -25);
}
TEST_CASE("Raw exponent digit quotas are checked before exponent accumulation",
          "[verification][negative][IO4-01]") {
    using namespace uvlas;
    core::ResourcePlan equal_plan{64ULL * 1024 * 1024};
    io::SchemaDiagnostic diagnostic;
    const auto equal = scoped_payload(
        "1e0001", "Number", with_quota(io::SchemaQuota::ExactDigits, 4), equal_plan, &diagnostic);
    CHECK(equal.has_value());
    core::ResourcePlan short_plan{64ULL * 1024 * 1024};
    const auto short_result = scoped_payload(
        "1e00001", "Number", with_quota(io::SchemaQuota::ExactDigits, 4), short_plan, &diagnostic);
    REQUIRE_FALSE(short_result.has_value());
    CHECK(short_result.error().code == core::FailureCode::ResourceExceeded);
    const auto detail = diagnostic.detail();
    if (!detail) {
        FAIL("Exponent digit rejection needs fixed detail");
        return;
    }
    CHECK(detail->quota == io::SchemaQuota::ExactDigits);
    CHECK(detail->required == 5);
}
TEST_CASE("Every structural error is selected before exact physical value checks",
          "[verification][negative][IO4-01]") {
    using namespace uvlas;
    auto value = test::schema_adapter::root();
    value.field("model").field("materials").elements[0].field("young") =
        test::schema_adapter::parse("-1");
    value.field("output").set("unrecognized", test::schema_adapter::parse("0"));
    const auto validate = [&](const test::schema_adapter::JsonToken& candidate,
                              io::SchemaReason reason) {
        const auto wire = candidate.wire();
        const io::InputLimits syntax{
            1000000,
            256,
            100000,
            100000,
            1000000,
            1000000,
            1000000,
            {io::ParserAllocationMode::EstimateAdmission, std::numeric_limits<std::uint64_t>::max()}};
        core::ResourcePlan input_plan{64ULL * 1024ULL * 1024};
        std::optional<io::InputDiagnostic> input_diagnostic;
        auto input =
            io::read_input_document(std::as_bytes(std::span{wire.data(), wire.size()}),
                                    io::SourceId{23}, syntax, input_plan, input_diagnostic);
        if (!input) {
            FAIL_CHECK("Phase control must pass the source envelope");
            return;
        }
        core::ResourcePlan plan{64ULL * 1024ULL * 1024};
        io::SchemaDiagnostic diagnostic;
        const auto checked = io::validate_case_schema(*input, schema_limits(), plan, &diagnostic);
        CHECK_FALSE(checked.has_value());
        const auto detail = diagnostic.detail();
        if (!detail) {
            FAIL_CHECK("Phase failure needs detail");
            return;
        }
        CHECK(detail->reason == reason);
        CHECK(detail->actual_kind.has_value());
        CHECK(detail->expected_type_name().has_value());
    };
    validate(value, io::SchemaReason::UnknownField);
    std::ranges::reverse(value.members);
    validate(value, io::SchemaReason::UnknownField);
    value.field("output").erase("unrecognized");
    validate(value, io::SchemaReason::InvalidValue);
}
TEST_CASE("A private payload view cannot transfer qualification to another document",
          "[verification][negative][IO4-01]") {
    using namespace uvlas;
    constexpr std::string_view wire =
        R"({"schema":"uvlas.case/1","units":1,"model":null,"flow":null,"analysis":null,"observations":null,"numerics":null,"resources":null,"output":null})";
    const io::InputLimits syntax{
        1000000,
        256,
        100000,
        100000,
        1000000,
        1000000,
        1000000,
        {io::ParserAllocationMode::EstimateAdmission, std::numeric_limits<std::uint64_t>::max()}};
    core::ResourcePlan first_plan{64ULL * 1024ULL * 1024}, second_plan{64ULL * 1024ULL * 1024};
    std::optional<io::InputDiagnostic> source_diagnostic;
    auto first = io::read_input_document(std::as_bytes(std::span{wire.data(), wire.size()}),
                                         io::SourceId{1}, syntax, first_plan, source_diagnostic);
    auto second = io::read_input_document(std::as_bytes(std::span{wire.data(), wire.size()}),
                                          io::SourceId{1}, syntax, second_plan, source_diagnostic);
    if (!first || !second) {
        FAIL("Two independent envelope controls are required");
        return;
    }
    const auto view = second->view();
    if (!view) {
        FAIL("Source view must exist");
        return;
    }
    const auto object = view->root.object();
    if (!object) {
        FAIL("Source root must be an object");
        return;
    }
    const auto node = object->find("units");
    if (!node) {
        FAIL("Payload node must exist");
        return;
    }
    core::ResourcePlan plan{64ULL * 1024ULL * 1024};
    io::SchemaDiagnostic diagnostic;
    const auto cross_owner = io::detail::validate_schema_payload(
        *first, {*node, declared_type("Number"), 1, io::SchemaDimension{{0, 0, 0, 0}}, false},
        schema_limits(), plan, &diagnostic);
    REQUIRE_FALSE(cross_owner.has_value());
    CHECK(cross_owner.error().code == core::FailureCode::InvalidInput);
    const auto detail = diagnostic.detail();
    if (!detail) {
        FAIL("Cross-owner rejection must retain its cause");
        return;
    }
    CHECK(detail->reason == io::SchemaReason::InvalidOwner);
    const auto transferred = std::move(*first);
    const auto moved = io::validate_case_schema(*first, schema_limits(), plan, &diagnostic);
    CHECK_FALSE(moved.has_value());
    if (!moved)
        CHECK(moved.error().operation == core::Operation::PrepareModel);
    const auto moved_detail = diagnostic.detail();
    CHECK(moved_detail.has_value());
    if (moved_detail)
        CHECK(moved_detail->reason == io::SchemaReason::InvalidOwner);
    CHECK(transferred.view().has_value());
}
TEST_CASE("Private registry scopes retain complete diagnostic paths and precharge traversal",
          "[verification][negative][IO4-01]") {
    using namespace uvlas;
    auto fixture = test::schema_adapter::root();
    fixture.set("units", test::schema_adapter::parse("1"));
    fixture.field("model").field("materials").elements[0].field("young") =
        test::schema_adapter::parse("-1");
    const auto wire = fixture.wire();
    const io::InputLimits syntax{
        1000000,
        256,
        100000,
        100000,
        1000000,
        1000000,
        1000000,
        {io::ParserAllocationMode::EstimateAdmission, std::numeric_limits<std::uint64_t>::max()}};
    core::ResourcePlan input_plan{64ULL * 1024 * 1024};
    std::optional<io::InputDiagnostic> input_diagnostic;
    auto input = io::read_input_document(std::as_bytes(std::span{wire.data(), wire.size()}),
                                         io::SourceId{25}, syntax, input_plan, input_diagnostic);
    if (!input) {
        FAIL("Registry fixture must pass the source envelope");
        return;
    }
    const auto view = input->view();
    if (!view) {
        FAIL("Registry source view must exist");
        return;
    }
    const auto object = view->root.object();
    if (!object) {
        FAIL("Registry root must be an object");
        return;
    }
    const auto units = object->find("units");
    if (!units) {
        FAIL("Payload field must exist");
        return;
    }
    const io::detail::SchemaPayloadRequest request{*units, declared_type("Number"), 1,
                                                   io::SchemaDimension{{0, 0, 0, 0}}, true};
    core::ResourcePlan plan{64ULL * 1024 * 1024};
    io::SchemaDiagnostic diagnostic;
    const auto result =
        io::detail::validate_schema_payload(*input, request, schema_limits(), plan, &diagnostic);
    CHECK_FALSE(result.has_value());
    const auto detail = diagnostic.detail();
    if (!detail) {
        FAIL("Registry rejection must retain its cause");
        return;
    }
    CHECK(detail->reason == io::SchemaReason::InvalidValue);
    CHECK(detail->path_state == io::SchemaPathState::Complete);
    const auto path = diagnostic.path();
    CHECK(path.has_value());
    if (path)
        CHECK(*path == "/model/materials/0/young");
    core::ResourcePlan no_work_plan{64ULL * 1024 * 1024};
    const auto no_work = io::detail::validate_schema_payload(
        *input, request, with_quota(io::SchemaQuota::Work, 0), no_work_plan, &diagnostic);
    CHECK_FALSE(no_work.has_value());
    const auto work_detail = diagnostic.detail();
    if (!work_detail) {
        FAIL("Work rejection must retain its cause");
        return;
    }
    CHECK(work_detail->quota == io::SchemaQuota::Work);
    CHECK(work_detail->work == 0);
    CHECK(no_work_plan.allocated_bytes() == 0);
}
TEST_CASE("All five frozen root designs reach the typed product validator",
          "[verification][IO4-01]") {
    using namespace uvlas;
    for (const auto& example : test::case_schema_corpus::roots) {
        INFO(example.name);
        const io::InputLimits syntax{
            1000000,
            256,
            100000,
            100000,
            1000000,
            1000000,
            1000000,
            {io::ParserAllocationMode::EstimateAdmission, std::numeric_limits<std::uint64_t>::max()}};
        core::ResourcePlan source_plan{64ULL * 1024ULL * 1024};
        std::optional<io::InputDiagnostic> source_diagnostic;
        auto document = io::read_input_document(
            std::as_bytes(std::span{example.json.data(), example.json.size()}), io::SourceId{7},
            syntax, source_plan, source_diagnostic);
        if (!document) {
            FAIL_CHECK("Frozen root must pass syntax");
            continue;
        }
        core::ResourcePlan plan{64ULL * 1024ULL * 1024};
        io::SchemaDiagnostic diagnostic;
        auto checked = io::validate_case_schema(*document, schema_limits(), plan, &diagnostic);
        const auto detail = diagnostic.detail();
        if (detail) {
            INFO("reason " << static_cast<unsigned>(detail->reason));
            if (detail->node)
                INFO("node " << *detail->node);
            const auto path = diagnostic.path();
            if (path)
                INFO("path " << *path);
            const auto type = detail->expected_type_name();
            if (type)
                INFO("expected " << *type);
        }
        if (!checked) {
            FAIL_CHECK("Frozen root must pass the complete typed schema");
            continue;
        }
        const auto root = checked->root();
        CHECK(root.has_value());
        const auto identity = checked->source_identity();
        if (!identity) {
            FAIL_CHECK("Source identity must be owned");
            continue;
        }
        CHECK(identity->source == io::SourceId{7});
        const auto usage = checked->memory_report();
        if (!usage) {
            FAIL_CHECK("Accounting must be owned");
            continue;
        }
        CHECK(usage->owned_bytes > 0);
        CHECK(usage->work > 0);
    }
}

TEST_CASE("Fixture editing preserves decimal tokens instead of rounding them", "[unit][IO4-01]") {
    using namespace uvlas::test::schema_adapter;
    const auto source = R"({"a":1e-400,"b":1.000000000000000000000000000001,"c":-0,"d":1.0})";
    const auto edited = parse(source);
    CHECK(edited.wire() == source);
}

TEST_CASE("Complete physical states and all initial negative stages use the product entry",
          "[verification][negative][IO4-01]") {
    using namespace uvlas;
    const auto validate = [&](std::string_view context, std::string_view perturbation,
                              bool accepted) {
        const auto wire = test::schema_adapter::physical_case(context, perturbation).wire();
        const io::InputLimits syntax{
            1000000,
            256,
            100000,
            100000,
            1000000,
            1000000,
            1000000,
            {io::ParserAllocationMode::EstimateAdmission, std::numeric_limits<std::uint64_t>::max()}};
        core::ResourcePlan source_plan{64ULL * 1024ULL * 1024};
        std::optional<io::InputDiagnostic> source_diagnostic;
        auto document =
            io::read_input_document(std::as_bytes(std::span{wire.data(), wire.size()}),
                                    io::SourceId{19}, syntax, source_plan, source_diagnostic);
        if (!document) {
            FAIL_CHECK("The typed corpus wrapper must pass syntax");
            return;
        }
        core::ResourcePlan plan{64ULL * 1024ULL * 1024};
        io::SchemaDiagnostic diagnostic;
        const auto checked =
            io::validate_case_schema(*document, schema_limits(), plan, &diagnostic);
        if (checked.has_value() != accepted) {
            const auto detail = diagnostic.detail();
            if (detail) {
                INFO("reason " << static_cast<unsigned>(detail->reason));
                const auto path = diagnostic.path();
                if (path) {
                    INFO("path " << *path);
                    CHECK(checked.has_value() == accepted);
                } else {
                    CHECK(checked.has_value() == accepted);
                }
            } else {
                CHECK(checked.has_value() == accepted);
            }
        } else {
            CHECK(checked.has_value() == accepted);
        }
        if (!accepted && !checked)
            CHECK(checked.error().code == core::FailureCode::InvalidInput);
        CHECK(document->view().has_value());
    };
    for (const auto& example : test::case_schema_corpus::physical) {
        INFO(example.name);
        validate(example.context_json, example.perturbation_json, true);
    }
    for (const auto& example : test::case_schema_corpus::negatives) {
        INFO(example.name);
        // A prepare-only contradiction must remain accepted by the schema layer.
        validate(example.context_json, example.perturbation_json,
                 example.stage == test::case_schema_corpus::ExpectedStage::Prepare);
    }
}
