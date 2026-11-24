#include "visualization_series.hpp"

#include <catch2/catch_test_macros.hpp>

#include <charconv>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

using namespace uvlas;
using namespace uvlas::io::detail;
using Side = VisualizationLayerSide;
using Measure = VisualizationMeasure;
using Basis = VisualizationBasis;
using Contribution = VisualizationContribution;
using FailureCode = core::FailureCode;

constexpr std::uint64_t budget = 65536;
constexpr auto extent = std::numeric_limits<std::uint64_t>::max();
constexpr VisualizationLimits limits{3, 8, 24, 12, budget};

template <class T>
T take(std::optional<T> value) {
    if (!value)
        throw std::runtime_error("Expected a file");
    return std::move(*value);
}

template <class T>
T take(core::Result<T> value) {
    if (!value)
        throw std::runtime_error("Expected a series, code=" +
                                 std::to_string(static_cast<unsigned>(value.error().code)));
    return std::move(*value);
}

std::array<VisualizationPoint, 4> points(double layer_fraction) {
    return {{{{0, 0, 0}, {0, 0, 1}, {0, 0, 1}, {0, 1, 0}, 0.5, 0.25, layer_fraction},
             {{1, 0, 0}, {1, 0, 1}, {0, 0, 1}, {0, 1, 0}, 0.5, 0.25, layer_fraction},
             {{0, 1, 0}, {0, 0, 0}, {0, 0, 1}, {0, 1, 0}, 0.5, 0.25, layer_fraction},
             {{1, 1, 0}, {1, 0, 0}, {0, 0, 1}, {0, 1, 0}, 0.5, 0.25, layer_fraction}}};
}

struct Fixture {
    std::array<VisualizationPoint, 4> below_points = points(1);
    std::array<VisualizationPoint, 4> above_points = points(0);
    std::array<VisualizationTriangle, 2> triangles{{{0, 1, 2}, {1, 3, 2}}};
    std::array<VisualizationSymmetric, 4> below_stress{};
    std::array<VisualizationSymmetric, 4> above_stress{};
    std::array<VisualizationSymmetric, 4> below_strain{};
    std::array<VisualizationMatrix, 4> below_first{};
    std::array<VisualizationTensorField, 3> below_tensors{
        {{Measure::Cauchy, Basis::World, Contribution::Total,
          std::span<const VisualizationSymmetric>{below_stress}},
         {Measure::Green, Basis::Material, Contribution::Total,
          std::span<const VisualizationSymmetric>{below_strain}},
         {Measure::FirstPk, Basis::Material, Contribution::Elastic,
          std::span<const VisualizationMatrix>{below_first}}}};
    std::array<VisualizationTensorField, 1> above_tensors{
        {{Measure::Cauchy, Basis::World, Contribution::Total,
          std::span<const VisualizationSymmetric>{above_stress}}}};
    std::array<VisualizationPart, 2> parts{
        {{{core::ComponentId{extent}, core::LayerId{7}, Side::Upper},
          below_points,
          triangles,
          below_tensors},
         {{core::ComponentId{extent}, core::LayerId{8}, Side::Lower},
          above_points,
          triangles,
          above_tensors}}};
    std::array<VisualizationFrame, 3> frames{
        {{0.125, core::AcceptedStateId{extent - 2}, 1, parts},
         {std::nextafter(0.125, 1.0), core::AcceptedStateId{extent - 1}, 1, parts},
         {0.125, core::AcceptedStateId{extent}, extent, parts}}};

    Fixture() {
        below_stress.fill({1, 2, 3, 4, 5, 6});
        above_stress.fill({10, 20, 30, 40, 50, 60});
        below_strain.fill({0.25, 0.5, 0.125, 0, 0, 0});
        below_first.fill({1, 2, 3, 4, 5, 6, 7, 8, 9});
    }
};

std::string_view array_text(std::string_view xml, std::string_view name) {
    const auto token = "Name=\"" + std::string{name} + "\"";
    auto position = xml.find(token);
    REQUIRE(position != xml.npos);
    position = xml.find('>', position);
    REQUIRE(position != xml.npos);
    const auto close = xml.find("</DataArray>", position);
    REQUIRE(close != xml.npos);
    return xml.substr(position + 1, close - position - 1);
}

double scalar(std::string_view text) {
    double value = 0;
    const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), value);
    REQUIRE(error == std::errc{});
    REQUIRE(end == text.data() + text.size());
    return value;
}

TEST_CASE("VTU and PVD own the accepted sample series with distinct material sides",
          "[IO4-05][QA-06]") {
    std::uint64_t exact_usage = 0;
    const auto documents = [&] {
        const Fixture input;
        core::ResourcePlan plan{budget};
        auto result = take(encode_visualization_series(input.frames, limits, plan));
        exact_usage = plan.allocated_bytes();
        return result;
    }();
    REQUIRE(documents.size() == 8);
    CHECK_FALSE(documents.at(8));
    const auto first_pvd = take(documents.at(0));
    const auto below = take(documents.at(1));
    const auto above = take(documents.at(2));
    const auto next = take(documents.at(3));
    const auto resumed_pvd = take(documents.at(5));
    const auto resumed = take(documents.at(6));
    CHECK(first_pvd.name == "segment-1.pvd");
    CHECK(resumed_pvd.name == "segment-18446744073709551615.pvd");
    CHECK(below.name == "accepted-18446744073709551613-part-0.vtu");
    CHECK(above.name == "accepted-18446744073709551613-part-1.vtu");
    CHECK(resumed.name == "accepted-18446744073709551615-part-0.vtu");
    CHECK(first_pvd.text.contains("timestep=\"0.125\" group=\"\" part=\"0\" file=\"" +
                                  std::string{below.name} + "\""));
    CHECK(first_pvd.text.contains("timestep=\"0.125\" group=\"\" part=\"1\" file=\"" +
                                  std::string{above.name} + "\""));
    CHECK(first_pvd.text.contains(next.name));
    CHECK_FALSE(first_pvd.text.contains(resumed.name));
    CHECK(resumed_pvd.text.contains("timestep=\"0.125\" group=\"\" part=\"0\" file=\"" +
                                    std::string{resumed.name} + "\""));
    CHECK_FALSE(resumed_pvd.text.contains(below.name));

    CHECK(below.text.contains("<Piece NumberOfPoints=\"4\" NumberOfCells=\"2\">"));
    CHECK(array_text(below.text, "reference_position_m") == "0 0 0 1 0 0 0 1 0 1 1 0 ");
    CHECK(array_text(below.text, "current_position_m") == "0 0 1 1 0 1 0 0 0 1 0 0 ");
    CHECK(array_text(below.text, "current_position_m") ==
          array_text(above.text, "current_position_m"));
    CHECK(array_text(above.text, "reference_normal") == "0 0 1 0 0 1 0 0 1 0 0 1 ");
    CHECK(array_text(above.text, "current_normal") == "0 1 0 0 1 0 0 1 0 0 1 0 ");
    CHECK(array_text(below.text, "reference_thickness_m") == "0.5 0.5 0.5 0.5 ");
    CHECK(array_text(below.text, "current_thickness_m") == "0.25 0.25 0.25 0.25 ");
    CHECK(array_text(below.text, "layer_fraction") == "1 1 1 1 ");
    CHECK(array_text(above.text, "layer_fraction") == "0 0 0 0 ");
    CHECK(array_text(below.text, "layer_side") == "2");
    CHECK(array_text(above.text, "layer_side") == "1");
    CHECK(array_text(below.text, "layer_id") == "7");
    CHECK(array_text(above.text, "layer_id") == "8");
    CHECK(array_text(below.text, "accepted_state_id") == "18446744073709551613");
    CHECK(array_text(resumed.text, "segment") == "18446744073709551615");
    CHECK(array_text(below.text, "component_id") == "18446744073709551615");
    CHECK(scalar(array_text(below.text, "time_s")) == 0.125);
    CHECK(scalar(array_text(next.text, "time_s")) == std::nextafter(0.125, 1.0));
    CHECK(array_text(below.text, "connectivity") == "0 1 2 1 3 2 ");
    CHECK(array_text(below.text, "offsets") == "3 6 ");
    CHECK(array_text(below.text, "types") == "5 5 ");
    CHECK(array_text(below.text, "stress_cauchy_total_world_Pa").starts_with("1 6 5 6 2 4 5 4 3 "));
    CHECK(array_text(above.text, "stress_cauchy_total_world_Pa")
              .starts_with("10 60 50 60 20 40 50 40 30 "));
    CHECK(array_text(below.text, "strain_green_material_1")
              .starts_with("0.25 0 0 0 0.5 0 0 0 0.125 "));
    CHECK(array_text(below.text, "stress_first_pk_elastic_material_Pa")
              .starts_with("1 2 3 4 5 6 7 8 9 "));
    CHECK_FALSE(above.text.contains("strain_green_material_1"));
    CHECK_FALSE(above.text.contains("stress_first_pk_elastic_material_Pa"));

    const Fixture replay_input;
    core::ResourcePlan exact{exact_usage};
    const auto replay = take(encode_visualization_series(replay_input.frames, limits, exact));
    REQUIRE(replay.size() == documents.size());
    CHECK(exact.allocated_bytes() == exact_usage);
    std::size_t output_bytes = 0;
    for (std::size_t i = 0; i < replay.size(); ++i) {
        const auto expected = take(documents.at(i));
        const auto actual = take(replay.at(i));
        output_bytes += actual.text.size();
        CHECK(actual.name == expected.name);
        CHECK(actual.text == expected.text);
    }
    auto caps = limits;
    caps.maximum_output_bytes = output_bytes - 1;
    core::ResourcePlan less_xml{budget};
    const auto too_many_bytes = encode_visualization_series(replay_input.frames, caps, less_xml);
    REQUIRE_FALSE(too_many_bytes);
    CHECK(too_many_bytes.error().code == FailureCode::ResourceExceeded);
    CHECK(less_xml.allocated_bytes() == 0);
    core::ResourcePlan less_total{exact_usage - 1};
    const auto too_large = encode_visualization_series(replay_input.frames, limits, less_total);
    REQUIRE_FALSE(too_large);
    CHECK(too_large.error().code == FailureCode::ResourceExceeded);
    CHECK(less_total.allocated_bytes() == 0);
}

TEST_CASE("Invalid or unadmitted visualization inputs produce no partial files",
          "[IO4-05][negative]") {
    Fixture input;
    auto caps = limits;
    auto expected = FailureCode::InvalidInput;
    // NOLINTNEXTLINE(misc-const-correctness): MSVC request_stop() requires a mutable source.
    std::stop_source stop;
    SECTION("time may repeat between segments only") {
        input.frames[1].time_seconds = input.frames[0].time_seconds;
    }
    SECTION("accepted identity must increase globally") {
        input.frames[2].accepted = input.frames[1].accepted;
    }
    SECTION("segments cannot be revisited") {
        input.frames[0].segment = 2;
    }
    SECTION("parts keep their identity within the same segment") {
        auto alien = input.parts;
        alien[1].key.layer = core::LayerId{9};
        input.frames[1].parts = alien;
        core::ResourcePlan plan{budget};
        const auto rejected = encode_visualization_series(input.frames, caps, plan);
        REQUIRE_FALSE(rejected);
        CHECK(rejected.error().code == expected);
        CHECK(plan.allocated_bytes() == 0);
        return;
    }
    SECTION("entries for the same side cannot overwrite each other") {
        input.parts[1].key = input.parts[0].key;
        input.parts[1].points = input.parts[0].points;
    }
    SECTION("a lower side cannot carry the upper endpoint") {
        input.above_points[0].layer_fraction = 1;
    }
    SECTION("thickness must be positive") {
        input.below_points[0].current_thickness_metres = 0;
    }
    SECTION("normal is part of the required geometry") {
        input.below_points[0].current_normal = {};
    }
    SECTION("point index is local to a piece") {
        input.triangles[0][2] = 4;
    }
    SECTION("collapsed triangle indices") {
        input.triangles[0][2] = input.triangles[0][1];
    }
    SECTION("first PK requires both off-diagonal terms") {
        input.below_tensors[2].values = std::span<const VisualizationSymmetric>{input.below_stress};
    }
    SECTION("tensor sample count must match geometry") {
        input.below_tensors[2].values =
            std::span<const VisualizationMatrix>{input.below_first}.first(3);
    }
    SECTION("Green strain does not have a viscous contribution") {
        input.below_tensors[1].contribution = Contribution::Viscous;
    }
    SECTION("non-finite coordinate") {
        input.below_points[0].current_metres[0] = std::numeric_limits<double>::infinity();
        expected = FailureCode::Nonfinite;
    }
    SECTION("non-finite tensor") {
        input.above_stress[3][5] = std::numeric_limits<double>::quiet_NaN();
        expected = FailureCode::Nonfinite;
    }
    SECTION("frame cap") {
        --caps.maximum_frames;
        expected = FailureCode::ResourceExceeded;
    }
    SECTION("file cap includes the PVD files") {
        --caps.maximum_files;
        expected = FailureCode::ResourceExceeded;
    }
    SECTION("point cap includes both interface sides") {
        --caps.maximum_points;
        expected = FailureCode::ResourceExceeded;
    }
    SECTION("triangle cap is cumulative") {
        --caps.maximum_triangles;
        expected = FailureCode::ResourceExceeded;
    }
    SECTION("cancelled") {
        REQUIRE(stop.request_stop());
        expected = FailureCode::Cancelled;
    }
    core::ResourcePlan plan{budget};
    const auto rejected = encode_visualization_series(input.frames, caps, plan, stop.get_token());
    REQUIRE_FALSE(rejected);
    CHECK(rejected.error().code == expected);
    CHECK(rejected.error().operation == core::Operation::WriteResults);
    CHECK(plan.allocated_bytes() == 0);
}

TEST_CASE("Empty visualization input owns no fabricated frames", "[IO4-05]") {
    core::ResourcePlan plan{0};
    const auto empty = take(encode_visualization_series({}, {0, 0, 0, 0, 0}, plan));
    CHECK(empty.size() == 0);
    CHECK_FALSE(empty.at(0));
    CHECK(plan.allocated_bytes() == 0);
}

} // namespace
