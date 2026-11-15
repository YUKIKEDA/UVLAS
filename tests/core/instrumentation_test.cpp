#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

#include <array>
#include <bit>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <uvlas/core/instrumentation.hpp>

namespace {
using namespace uvlas::core;
constexpr Operation operation = Operation::PrepareModel;
} // namespace

TEST_CASE("Phase samples and disabled numerical bytes", "[unit][N4-09]") {
    StageRecorder recorder{true, operation};
    REQUIRE(recorder.record_sample(MeasurementStage::Structure, 13));
    REQUIRE(recorder.record_sample(MeasurementStage::Structure, 29));
    const auto times = recorder.snapshot();
    REQUIRE(times);
    CHECK((*times)[1].calls == 2);
    CHECK((*times)[1].elapsed_nanoseconds == 42);
    for (const bool enabled : {false, true}) {
        StageRecorder scoped{enabled, operation};
        std::array<double, 4> values{1.0, -2.0, 0.0, 0.5};
        {
            const StageScope scope{scoped, MeasurementStage::LinearSolve};
            for (double& value : values) {
                value *= 2;
            }
        }
        const std::array<double, 4> expected{2.0, -4.0, 0.0, 1.0};
        CHECK(std::bit_cast<std::array<std::byte, sizeof(values)>>(values) ==
              std::bit_cast<std::array<std::byte, sizeof(expected)>>(expected));
        const auto result = scoped.snapshot();
        REQUIRE(result);
        CHECK((*result)[3].calls == (enabled ? 1 : 0));
        if (!enabled) {
            CHECK((*result)[3].elapsed_nanoseconds == 0);
        }
    }
}

TEST_CASE("All measurement phases preserve exact samples", "[unit][N4-09]") {
    StageRecorder recorder{true, operation};
    constexpr std::array<std::pair<MeasurementStage, std::string_view>, 8> phases{
        {{MeasurementStage::Preparation, "preparation"},
         {MeasurementStage::Structure, "structure"},
         {MeasurementStage::Aerodynamics, "aerodynamics"},
         {MeasurementStage::LinearSolve, "linear_solve"},
         {MeasurementStage::NonlinearSolve, "nonlinear_solve"},
         {MeasurementStage::TimeStep, "time_step"},
         {MeasurementStage::FlutterSearch, "flutter_search"},
         {MeasurementStage::Checkpoint, "checkpoint"}}};
    constexpr auto count = phases.size();
    for (std::size_t i = 0; i < count; ++i) {
        REQUIRE(recorder.record_sample(phases[i].first, i + 1));
    }
    const auto times = recorder.snapshot();
    REQUIRE(times);
    const ResourcePlan plan{0};
    const MeasurementConditions conditions{{"fixture", "fixture-cpu", 1, 1}, "all-stages", 1, "{}"};
    const auto encoded =
        measurement_report_json(conditions, recorder, plan, ProcessMemory{1, 1}, operation);
    REQUIRE(encoded);
    const auto report = nlohmann::json::parse(*encoded);
    for (std::size_t i = 0; i < count; ++i) {
        const auto index = static_cast<std::size_t>(phases[i].first);
        CHECK((*times)[index].calls == 1);
        CHECK((*times)[index].elapsed_nanoseconds == i + 1);
        CHECK(std::string_view{
                  report.at("stages").at(index).at("stage").get_ref<const std::string&>()} ==
              phases[i].second);
        CHECK(report.at("stages").at(index).at("calls") == 1);
        CHECK(report.at("stages").at(index).at("elapsed_nanoseconds") == i + 1);
    }
}

TEST_CASE("Phase unwinding and retained failures", "[negative][N4-09]") {
    StageRecorder recorder{true, operation};
    bool caught = false;
    try {
        const StageScope scope{recorder, MeasurementStage::Structure};
        throw std::runtime_error{"cancel"};
    } catch (const std::runtime_error&) {
        caught = true;
    }
    CHECK(caught);
    const auto times = recorder.snapshot();
    REQUIRE(times);
    CHECK((*times)[1].calls == 1);
    REQUIRE(recorder.record_sample(MeasurementStage::Preparation,
                                   std::numeric_limits<std::uint64_t>::max()));
    const auto overflow = recorder.record_sample(MeasurementStage::Preparation, 1);
    REQUIRE_FALSE(overflow);
    CHECK(overflow.error().code == FailureCode::ResourceExceeded);
    CHECK_FALSE(recorder.snapshot());
    StageRecorder invalid{true, operation};
    const auto bad_stage = invalid.record_sample(MeasurementStage::Count, 0);
    REQUIRE_FALSE(bad_stage);
    CHECK(bad_stage.error().code == FailureCode::InvalidInput);
    CHECK_FALSE(invalid.snapshot());
}

TEST_CASE("Measurement JSON separates RSS and payload", "[integration][N4-09]") {
    StageRecorder recorder{true, operation};
    REQUIRE(recorder.record_sample(MeasurementStage::Structure, 42));
    ResourcePlan resources{64};
    const auto buffer = resources.allocate<double>(2, operation);
    REQUIRE(buffer);
    const auto memory = sample_process_memory(operation);
    REQUIRE(memory);
    CHECK(memory->resident_bytes > 0);
    CHECK(memory->peak_resident_bytes > 0);
    const MachineDescription machine{"test-machine", "test-cpu", 8, 16ULL * 1024 * 1024 * 1024};
    const MeasurementConditions conditions{machine, "phase-fixture", 2, R"({"n":2})"};
    const auto encoded =
        measurement_report_json(conditions, recorder, resources, *memory, operation);
    REQUIRE(encoded);
    const auto report = nlohmann::json::parse(*encoded);
    CHECK(report.at("schema") == "uvlas.measurement/1");
    CHECK(report.at("threads") == 2);
    CHECK(report.at("settings").at("n") == 2);
    CHECK(report.at("memory").at("planned_payload_bytes") == 16);
    CHECK(report.at("memory").at("peak_resident_bytes") == memory->peak_resident_bytes);
    CHECK(report.at("build").contains("instruction_set"));
    CHECK(report.at("build").contains("commit"));
    CHECK(report.at("stages").at(1).at("elapsed_nanoseconds") == 42);
    CHECK(report.at("machine").at("cpu") == "test-cpu");
    CHECK(report.at("machine").at("logical_processors") == 8);
    for (const auto& bad : {MeasurementConditions{{"", "cpu", 1, 1}, "case", 1, "{}"},
                            MeasurementConditions{{"machine", "", 1, 1}, "case", 1, "{}"},
                            MeasurementConditions{{"machine", "cpu", 0, 1}, "case", 1, "{}"},
                            MeasurementConditions{{"machine", "cpu", 1, 0}, "case", 1, "{}"},
                            MeasurementConditions{machine, "", 1, "{}"},
                            MeasurementConditions{machine, "case", 0, "{}"},
                            MeasurementConditions{machine, "case", 1, "[]"},
                            MeasurementConditions{machine, "case", 1, "{"}}) {
        const auto rejected = measurement_report_json(bad, recorder, resources, *memory, operation);
        REQUIRE_FALSE(rejected);
        CHECK(rejected.error().code == FailureCode::InvalidInput);
    }
}
