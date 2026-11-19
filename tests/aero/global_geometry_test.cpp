#include "aero/global_geometry.hpp"
#include "geometry_snapshot_fixture.hpp"
#include "global_geometry_fixture.hpp"
#include "global_kinematics_fixture.hpp"

#include <catch2/catch_test_macros.hpp>

#include <fstream>
#include <future>
#include <limits>

namespace {
using namespace uvlas;
namespace s = aero::detail;
constexpr std::uint64_t budget = 64ULL * 1024 * 1024;
s::GlobalLimits limits() {
    return {std::numeric_limits<std::uint64_t>::max() / 4,
            numerics::linear::EstimateAdmission{4096}};
}
template <class T>
T take(s::GlobalResult<T> value) {
    if (!value) {
        FAIL_CHECK("global failure " << static_cast<unsigned>(value.error().origin.cause.code));
        throw std::runtime_error("global operation failed");
    }
    return std::move(*value);
}
std::unique_ptr<s::GlobalWorkspace> workspace(core::ResourcePlan& plan) {
    auto w = s::GlobalWorkspace::create(plan, numerics::linear::EstimateAdmission{4096});
    REQUIRE(w);
    return std::move(*w);
}
nlohmann::json reference_case(std::string_view name) {
    std::ifstream stream{UVLAS_GLOBAL_GEOMETRY_FIXTURE};
    nlohmann::json data;
    stream >> data;
    for (const auto& row : data.at("cases"))
        if (row.at("name").get<std::string>() == name)
            return row;
    throw std::runtime_error("missing fixture");
}
void contact_budgets(const s::GlobalTopology& topology, const s::OwnedGeometry& geometry,
                     s::GlobalWorkspace& work, std::uint64_t facets,
                     std::uint64_t expected_candidates) {
    core::ResourcePlan measured{budget};
    auto state = take(s::check_global_geometry(topology, geometry, limits(), work, measured));
    const auto planned = take(state.work_plan());
    const auto nodes = 2 * facets - 1, pairs = nodes * (nodes + 1) / 2;
    constexpr std::uint64_t limb_count = 1024, decode = 8 * limb_count, comparison = 3 * limb_count;
    const auto expected_traversal = 2 * pairs * (6 * comparison + 24 * decode + 3 + 1 + 1);
    CHECK(planned.candidates == expected_candidates);
    CHECK(planned.traversal_reserved_steps == expected_traversal);
    CHECK(planned.traversal_executed_steps > 0);
    CHECK(planned.traversal_executed_steps <= expected_traversal);
    CHECK(planned.traversal_executed_steps % 2 == 0);
    auto exact = limits();
    exact.max_steps = planned.reserved_steps;
    core::ResourcePlan equal{planned.payload_bytes};
    CHECK(s::check_global_geometry(topology, geometry, exact, work, equal));
    --exact.max_steps;
    core::ResourcePlan short_work{budget};
    auto rejected = s::check_global_geometry(topology, geometry, exact, work, short_work);
    REQUIRE_FALSE(rejected);
    CHECK(rejected.error().origin.cause.code == core::FailureCode::ResourceExceeded);
    const auto resource = rejected.error().origin.cause.details.resources;
    if (!resource) {
        FAIL("Missing work admission diagnostic");
        return;
    }
    CHECK(resource->kind == core::ResourceKind::Iterations);
    CHECK(resource->required == planned.reserved_steps);
    core::ResourcePlan short_bytes{planned.payload_bytes - 1};
    CHECK_FALSE(s::check_global_geometry(topology, geometry, limits(), work, short_bytes));
    CHECK(short_bytes.allocated_bytes() == 0);
}
} // namespace
TEST_CASE("Physical aliases compare all columns time terms and exact velocity assertions",
          "[UW-D07][UI-D03][N4-02][N4-03][N4-06]") {
    test_global_geometry::KinematicsFixture fixture;
    core::ResourcePlan plan{budget};
    auto old = test_geometry::workspace(plan);
    auto w = workspace(plan);
    auto geometry = test_geometry::take(
        s::prepare_owned_geometry(fixture.geometry(), test_geometry::limits(), *old, plan));
    auto baseline = take(s::check_alias_kinematics(geometry, fixture.input(), limits(), *w, plan));
    const auto owner = take(baseline.velocity(fixture.base.points[1], 0)),
               alias = take(baseline.velocity(s::PointId{0}, 0));
    CHECK(owner.limbs.data() == alias.limbs.data());
    CHECK(owner.sign == alias.sign);
    CHECK(owner.first_limb == alias.first_limb);
    for (std::size_t column = 0; column < 3; ++column) {
        fixture.alias_coefficients[column].value = s::DoubleWord{1, 0x1p-100};
        auto result = s::check_alias_kinematics(geometry, fixture.input(), limits(), *w, plan);
        REQUIRE_FALSE(result);
        CHECK(result.error().origin.cause.code == core::FailureCode::InvalidInput);
        CHECK(result.error().origin.problem == s::GeometryProblem{s::PointId{0}});
        fixture.alias_coefficients[column].value = s::SingleWord{1};
    }
    // Even v=0 must not turn a full-J mismatch into a qualified one-direction JVP.
    for (auto& dof : fixture.dofs)
        dof.velocity = s::SingleWord{0};
    fixture.alias_coefficients.back().value = s::SingleWord{2};
    CHECK_FALSE(s::check_alias_kinematics(geometry, fixture.input(), limits(), *w, plan));
    fixture.alias_coefficients.back().value = s::SingleWord{1};
    for (std::size_t xyz = 0; xyz < 3; ++xyz) {
        fixture.points.back().explicit_time[xyz] = s::SingleWord{0x1p-100};
        CHECK_FALSE(s::check_alias_kinematics(geometry, fixture.input(), limits(), *w, plan));
        fixture.points.back().explicit_time[xyz] = s::SingleWord{0};
    }
    fixture.dofs[0].velocity = s::SingleWord{1};
    fixture.dofs[1].velocity = s::SingleWord{0x1p-200};
    fixture.dofs[2].velocity = s::SingleWord{0x1p-400};
    std::array<s::ExactVelocityInput, 3> assertions{s::ExactVelocityInput{0, 0, {}},
                                                    s::ExactVelocityInput{0, 0, {}},
                                                    s::ExactVelocityInput{0, 0, {}}};
    for (std::size_t k = 0; k < 3; ++k) {
        const auto exact = take(baseline.velocity(s::PointId{0}, k));
        assertions[k] = {exact.sign, exact.first_limb, exact.limbs};
    }
    fixture.points.back().velocity = assertions;
    CHECK(s::check_alias_kinematics(geometry, fixture.input(), limits(), *w, plan));
    const std::array<std::uint32_t, 1> wrong{1};
    for (std::size_t k = 0; k < 3; ++k) {
        auto inconsistent = assertions;
        inconsistent[k] = {1, 0, wrong};
        fixture.points.back().velocity = inconsistent;
        CHECK_FALSE(s::check_alias_kinematics(geometry, fixture.input(), limits(), *w, plan));
    }
    fixture.points.back().velocity.reset();
    std::reverse(fixture.points.begin(), fixture.points.end());
    CHECK(s::check_alias_kinematics(geometry, fixture.input(), limits(), *w, plan));
    fixture.dofs[0].velocity = s::SingleWord{std::numeric_limits<double>::infinity()};
    auto nonfinite = s::check_alias_kinematics(geometry, fixture.input(), limits(), *w, plan);
    REQUIRE_FALSE(nonfinite);
    CHECK(nonfinite.error().origin.cause.code == core::FailureCode::Nonfinite);
    CHECK(take(baseline.velocity(s::PointId{0}, 0)).limbs[0] == owner.limbs[0]);
}
TEST_CASE("Global geometry rejects zero area but admits a proper half-turn without a provider "
          "certificate",
          "[UW-D07][UI-D03][N4-03]") {
    test_geometry::TriangleFixture fixture;
    fixture.current[7] = -fixture.reference[7];
    core::ResourcePlan plan{budget};
    auto old = test_geometry::workspace(plan);
    auto work = workspace(plan);
    auto rotated = test_geometry::take(
        s::prepare_owned_geometry(fixture.input(), test_geometry::limits(), *old, plan));
    const std::array roles{s::FaceRoleInput{fixture.faces[0].id, s::AeroFaceRole::Bound}};
    const std::array edges{s::EdgeBindingInput{fixture.edges[0].id, s::CanonicalEdge{}},
                           s::EdgeBindingInput{fixture.edges[1].id, s::CanonicalEdge{}},
                           s::EdgeBindingInput{fixture.edges[2].id, s::CanonicalEdge{}}};
    auto topology =
        take(s::build_global_topology(rotated, {roles, edges, {}}, limits(), *work, plan));
    CHECK(s::check_global_geometry(topology, rotated, limits(), *work, plan));
    static_assert(std::is_same_v<decltype(s::GlobalContactState::mapping_status()),
                                 s::PhysicalMappingNotChecked>);
    fixture.current[6] = 6;
    fixture.current[7] = 0;
    fixture.stamp.current = core::GeometryGeneration{44};
    // A1 already rejects a zero-area source before any global qualification can
    // be requested; reuse that mandatory predecessor rather than forge an owner.
    auto collapsed =
        s::prepare_owned_geometry(fixture.input(), test_geometry::limits(), *old, plan);
    REQUIRE_FALSE(collapsed);
    CHECK(collapsed.error().cause.code == core::FailureCode::GeometryInvalid);
    CHECK(collapsed.error().problem == s::GeometryProblem{s::FacetId{70}});
}
TEST_CASE("Sparse alias rows reject duplicate unsorted unknown and explicit-zero columns",
          "[UW-D07][UI-D03][N4-02]") {
    test_global_geometry::KinematicsFixture fixture;
    core::ResourcePlan plan{budget};
    auto old = test_geometry::workspace(plan);
    auto work = workspace(plan);
    auto geometry = test_geometry::take(
        s::prepare_owned_geometry(fixture.geometry(), test_geometry::limits(), *old, plan));
    fixture.alias_coefficients[1].column = core::DofId{1};
    CHECK_FALSE(s::check_alias_kinematics(geometry, fixture.input(), limits(), *work, plan));
    fixture.alias_coefficients[1].column = core::DofId{0};
    CHECK_FALSE(s::check_alias_kinematics(geometry, fixture.input(), limits(), *work, plan));
    fixture.alias_coefficients[1].column = core::DofId{2};
    fixture.alias_coefficients[2].column = core::DofId{999};
    auto unknown = s::check_alias_kinematics(geometry, fixture.input(), limits(), *work, plan);
    REQUIRE_FALSE(unknown);
    CHECK(unknown.error().origin.cause.details.problem == core::ProblemId{core::DofId{999}});
    fixture.alias_coefficients[2].column = core::DofId{3};
    fixture.alias_coefficients[2].value = s::SingleWord{0};
    fixture.coefficients[2].value = s::SingleWord{0};
    CHECK_FALSE(s::check_alias_kinematics(geometry, fixture.input(), limits(), *work, plan));
}
TEST_CASE("Reversed edge aliases preserve canonical incidence and do not legalize internal "
          "trailing seams",
          "[UW-D04][UW-D07][UW-D08][UI-D03]") {
    test_global_geometry::PairFixture fixture{reference_case("coplanar-edge")};
    // Raw alias IDs ascend while their physical owners descend.
    for (const auto& link :
         std::array<std::pair<std::uint64_t, std::uint64_t>, 2>{{{20, 2}, {30, 1}}}) {
        const auto found =
            std::find(fixture.points.begin(), fixture.points.end(), s::PointId{link.second});
        REQUIRE(found != fixture.points.end());
        const auto ordinal = static_cast<std::size_t>(found - fixture.points.begin());
        for (std::size_t k = 0; k < 3; ++k) {
            fixture.high.push_back(fixture.high[3 * ordinal + k]);
            fixture.low.push_back(fixture.low[3 * ordinal + k]);
        }
        fixture.points.emplace_back(link.first);
        fixture.bindings.emplace_back(s::ReferenceMapId{1}, s::MappingGeneration{1},
                                      s::AliasOf{s::PointId{link.second}});
    }
    fixture.triangles[1].points[0] = s::PointId{20};
    fixture.triangles[1].points[1] = s::PointId{30};
    fixture.edges.emplace_back(s::EdgeId{200}, std::array{s::PointId{20}, s::PointId{30}});
    fixture.edge_bindings.push_back({s::EdgeId{200}, s::EdgeId{10}});
    for (auto& edge : fixture.triangles[1].edges)
        if (edge == s::EdgeId{10})
            edge = s::EdgeId{200};
    core::ResourcePlan plan{budget};
    auto old = test_geometry::workspace(plan);
    auto w = workspace(plan);
    auto geometry = test_geometry::take(
        s::prepare_owned_geometry(fixture.input(), test_geometry::limits(), *old, plan));
    const std::array te{s::TrailingEdgeInput{s::EdgeId{200}, {s::PointId{20}, s::PointId{30}}}};
    const auto rejected = s::build_global_topology(
        geometry, {fixture.roles, fixture.edge_bindings, te}, limits(), *w, plan);
    REQUIRE_FALSE(rejected);
    CHECK(rejected.error().origin.cause.code == core::FailureCode::InvalidInput);
    fixture.roles[1].role = s::AeroFaceRole::Wake;
    auto topology = take(s::build_global_topology(
        geometry, {fixture.roles, fixture.edge_bindings, te}, limits(), *w, plan));
    REQUIRE(topology.trailing_edges().size() == 1);
    CHECK(topology.trailing_edges()[0].edge == s::EdgeId{10});
    CHECK(topology.trailing_edges()[0].coefficient == -1);
    CHECK(take(topology.release_key(s::PointId{20}, s::ReleaseEpoch{7})) ==
          take(topology.release_key(s::PointId{2}, s::ReleaseEpoch{7})));
    CHECK(take(s::check_global_geometry(topology, geometry, limits(), *w, plan)).active());
    fixture.edge_bindings.back().ownership = s::EdgeId{200};
    CHECK_FALSE(s::build_global_topology(geometry, {fixture.roles, fixture.edge_bindings, te},
                                         limits(), *w, plan));
}
TEST_CASE(
    "Global contact results preserve exact input content and are independent across thread counts",
    "[UW-D07][N4-03][N4-08][N4-09]") {
    test_global_geometry::PairFixture fixture{reference_case("coplanar-edge")};
    core::ResourcePlan plan{budget};
    auto old = test_geometry::workspace(plan);
    auto w = workspace(plan);
    auto geometry = test_geometry::take(
        s::prepare_owned_geometry(fixture.input(), test_geometry::limits(), *old, plan));
    auto topology =
        take(s::build_global_topology(geometry, fixture.topology(), limits(), *w, plan));
    const auto original = take(s::check_global_geometry(topology, geometry, limits(), *w, plan));
    const auto expected = take(original.work_plan());
    for (const auto workers : std::array<std::size_t, 3>{1, 2, 4}) {
        std::vector<std::future<s::GlobalWorkPlan>> tasks;
        tasks.reserve(workers);
        for (std::size_t i = 0; i < workers; ++i)
            tasks.push_back(std::async(std::launch::async, [&] {
                core::ResourcePlan memory{budget};
                auto own =
                    s::GlobalWorkspace::create(memory, numerics::linear::EstimateAdmission{4096});
                if (!own)
                    throw std::runtime_error("workspace failed");
                auto result = s::check_global_geometry(topology, geometry, limits(), **own, memory);
                if (!result)
                    throw std::runtime_error("contact failed");
                return *result->work_plan();
            }));
        for (auto& task : tasks) {
            const auto actual = task.get();
            CHECK(actual.payload_bytes == expected.payload_bytes);
            CHECK(actual.reserved_steps == expected.reserved_steps);
            CHECK(actual.executed_steps == expected.executed_steps);
            CHECK(actual.candidates == expected.candidates);
            CHECK(actual.metadata_forecast == expected.metadata_forecast);
        }
    }
    const auto identity = test_geometry::take(geometry.identity());
    fixture.high[0] += 0.125;
    auto changed = test_geometry::take(
        s::prepare_owned_geometry(fixture.input(), test_geometry::limits(), *old, plan));
    auto wrong = s::check_global_geometry(topology, changed, limits(), *w, plan);
    REQUIRE_FALSE(wrong);
    CHECK(wrong.error().origin.cause.code == core::FailureCode::InvalidInput);
    // MSVC request_stop() is non-const; libstdc++ declares it const.
    // NOLINTNEXTLINE(misc-const-correctness)
    std::stop_source stop;
    stop.request_stop();
    core::ResourcePlan unused{budget};
    auto cancelled =
        s::check_global_geometry(topology, geometry, limits(), *w, unused, stop.get_token());
    REQUIRE_FALSE(cancelled);
    CHECK(cancelled.error().origin.cause.code == core::FailureCode::Cancelled);
    CHECK(unused.allocated_bytes() == 0);
    CHECK(test_geometry::take(geometry.identity()) == identity);
}
TEST_CASE("Global contacts match the independent feasible-polytope literals under all permutations",
          "[UW-D07][UI-D03][N4-02][N4-06]") {
    std::ifstream stream{UVLAS_GLOBAL_GEOMETRY_FIXTURE};
    REQUIRE(stream.good());
    nlohmann::json reference;
    stream >> reference;
    core::ResourcePlan work_memory{budget};
    auto old = test_geometry::workspace(work_memory);
    auto w = workspace(work_memory);
    for (const auto& record : reference.at("cases")) {
        for (const auto& setting : std::array<std::pair<int, bool>, 4>{
                 {{0, false}, {-900, false}, {900, false}, {0, true}}}) {
            std::array<std::size_t, 3> pa{0, 1, 2};
            do {
                std::array<std::size_t, 3> pb{0, 1, 2};
                do {
                    INFO(record.at("name").get<std::string>()
                         << " permutations " << pa[0] << pa[1] << pa[2] << " " << pb[0] << pb[1]
                         << pb[2] << " scale " << setting.first << " origin " << setting.second);
                    test_global_geometry::PairFixture const fixture{record, pa, pb, setting.first,
                                                                    setting.second};
                    core::ResourcePlan plan{budget};
                    auto g = test_geometry::take(s::prepare_owned_geometry(
                        fixture.input(), test_geometry::limits(), *old, plan));
                    auto topology =
                        take(s::build_global_topology(g, fixture.topology(), limits(), *w, plan));
                    auto checked = s::check_global_geometry(topology, g, limits(), *w, plan);
                    if (record.at("name") == "disjoint" &&
                        pa == std::array<std::size_t, 3>{0, 1, 2} && pb == pa &&
                        setting.first == 0 && !setting.second)
                        contact_budgets(topology, g, *w, 2, 0);
                    CHECK(checked.has_value() == record.at("accepted").get<bool>());
                    if (!checked) {
                        CHECK(checked.error().origin.cause.code ==
                              core::FailureCode::GeometryInvalid);
                        REQUIRE(checked.error().other_facet);
                        CHECK(checked.error().other_facet->value() == 71);
                    }
                } while (std::next_permutation(pb.begin(), pb.end()));
            } while (std::next_permutation(pa.begin(), pa.end()));
        }
    }
}
TEST_CASE("Global topology retains a single bound TE and its epoch owners",
          "[UW-D04][UW-D07][UW-D08][UI-D03][N4-03]") {
    test_geometry::TriangleFixture fixture;
    core::ResourcePlan plan{budget};
    auto old = test_geometry::workspace(plan);
    auto w = workspace(plan);
    auto geometry = test_geometry::take(
        s::prepare_owned_geometry(fixture.input(), test_geometry::limits(), *old, plan));
    const std::array roles{s::FaceRoleInput{fixture.faces[0].id, s::AeroFaceRole::Bound}};
    const std::array edges{s::EdgeBindingInput{fixture.edges[0].id, s::CanonicalEdge{}},
                           s::EdgeBindingInput{fixture.edges[1].id, s::CanonicalEdge{}},
                           s::EdgeBindingInput{fixture.edges[2].id, s::CanonicalEdge{}}};
    const std::array te{
        s::TrailingEdgeInput{fixture.edges[0].id, {fixture.points[1], fixture.points[0]}}};
    auto topology =
        take(s::build_global_topology(geometry, {roles, edges, te}, limits(), *w, plan));
    REQUIRE(topology.trailing_edges().size() == 1);
    CHECK(topology.trailing_edges()[0].coefficient == -1);
    CHECK(topology.trailing_point_count() == 2);
    const auto a = take(topology.release_key(fixture.points[0], s::ReleaseEpoch{0}));
    const auto b = take(topology.release_key(fixture.points[0], s::ReleaseEpoch{UINT64_MAX}));
    CHECK(a != b);
    auto checked = take(s::check_global_geometry(topology, geometry, limits(), *w, plan));
    CHECK(checked.active());
    CHECK(take(checked.work_plan()).candidates == 0);
    CHECK(take(topology.work_plan()).executed_steps <= take(topology.work_plan()).reserved_steps);
}
TEST_CASE("Alias kinematics retain three disjoint velocity bit bands",
          "[UW-D07][UI-D03][N4-02][N4-06]") {
    test_geometry::TriangleFixture fixture;
    core::ResourcePlan plan{budget};
    auto old = test_geometry::workspace(plan);
    auto w = workspace(plan);
    auto geometry = test_geometry::take(
        s::prepare_owned_geometry(fixture.input(), test_geometry::limits(), *old, plan));
    const std::array<s::ScalarWords, 3> axis{s::SingleWord{1}, s::SingleWord{0}, s::SingleWord{0}};
    const std::array dofs{
        s::GlobalDofInput{
            core::DofId{1}, {1, 0, 0, 0, 0, 0, 0}, s::DofBasisId{1}, axis, s::SingleWord{1}},
        s::GlobalDofInput{
            core::DofId{2}, {1, 0, 0, 0, 0, 0, 0}, s::DofBasisId{1}, axis, s::SingleWord{0x1p-200}},
        s::GlobalDofInput{core::DofId{3},
                          {1, 0, 0, 0, 0, 0, 0},
                          s::DofBasisId{1},
                          axis,
                          s::SingleWord{0x1p-400}}};
    const std::array entries{s::GlobalJacobianEntry{core::DofId{1}, s::SingleWord{1}},
                             s::GlobalJacobianEntry{core::DofId{2}, s::SingleWord{1}},
                             s::GlobalJacobianEntry{core::DofId{3}, s::SingleWord{1}}};
    const std::array<s::ScalarWords, 3> zero{s::SingleWord{0}, s::SingleWord{0}, s::SingleWord{0}};
    std::array points{s::PointKinematicsInput{fixture.points[0], {entries, {}, {}}, zero, {}},
                      s::PointKinematicsInput{fixture.points[1], {entries, {}, {}}, zero, {}},
                      s::PointKinematicsInput{fixture.points[2], {entries, {}, {}}, zero, {}}};
    auto state = take(s::check_alias_kinematics(geometry, {s::DofLayoutGeneration{0}, dofs, points},
                                                limits(), *w, plan));
    const auto v = take(state.velocity(fixture.points[0], 0));
    CHECK(v.sign == 1);
    CHECK(v.first_limb == 121);
    REQUIRE(v.limbs.size() == 14);
    CHECK(v.limbs[0] == (UINT32_C(1) << 24U));
    CHECK(v.limbs[7] == UINT32_C(1));
    CHECK(v.limbs[13] == (UINT32_C(1) << 8U));
}
TEST_CASE("Global qualifications outlive all input arrays workspaces and resource plans",
          "[UW-D07][N4-03][N4-08]") {
    const auto make_topology = []() {
        test_geometry::TriangleFixture fixture;
        core::ResourcePlan plan{budget};
        auto old = test_geometry::workspace(plan);
        auto work = workspace(plan);
        auto source = test_geometry::take(
            s::prepare_owned_geometry(fixture.input(), test_geometry::limits(), *old, plan));
        const std::array roles{s::FaceRoleInput{fixture.faces[0].id, s::AeroFaceRole::Bound}};
        const std::array edges{s::EdgeBindingInput{fixture.edges[0].id, s::CanonicalEdge{}},
                               s::EdgeBindingInput{fixture.edges[1].id, s::CanonicalEdge{}},
                               s::EdgeBindingInput{fixture.edges[2].id, s::CanonicalEdge{}}};
        const std::array te{
            s::TrailingEdgeInput{fixture.edges[0].id, {fixture.points[1], fixture.points[0]}}};
        return take(s::build_global_topology(source, {roles, edges, te}, limits(), *work, plan));
    };
    auto topology = make_topology();
    CHECK(topology.active());
    CHECK(take(topology.geometry()).get().active());
    CHECK(take(topology.release_key(s::PointId{11}, s::ReleaseEpoch{4})).point == s::PointId{11});
    auto moved = std::move(topology);
    // NOLINTBEGIN(bugprone-use-after-move): testing the documented inactive moved-from state.
    CHECK_FALSE(topology.active());
    CHECK_FALSE(topology.geometry());
    // NOLINTEND(bugprone-use-after-move)
    CHECK(moved.active());
    auto contact = [&]() {
        core::ResourcePlan plan{budget};
        auto work = workspace(plan);
        return take(
            s::check_global_geometry(moved, take(moved.geometry()).get(), limits(), *work, plan));
    }();
    CHECK(take(contact.geometry()).get().active());
    auto kinematics = []() {
        test_global_geometry::KinematicsFixture const fixture;
        core::ResourcePlan plan{budget};
        auto old = test_geometry::workspace(plan);
        auto work = workspace(plan);
        auto source = test_geometry::take(
            s::prepare_owned_geometry(fixture.geometry(), test_geometry::limits(), *old, plan));
        return take(s::check_alias_kinematics(source, fixture.input(), limits(), *work, plan));
    }();
    CHECK(take(kinematics.velocity(s::PointId{0}, 0)).limbs.size() == 14);
    CHECK(take(kinematics.geometry()).get().active());
    auto moved_kinematics = std::move(kinematics);
    // NOLINTBEGIN(bugprone-use-after-move): moved-from access must fail by value.
    CHECK_FALSE(kinematics.active());
    CHECK_FALSE(kinematics.velocity(s::PointId{0}, 0));
    // NOLINTEND(bugprone-use-after-move)
    CHECK(moved_kinematics.active());
    auto moved_contact = std::move(contact);
    // NOLINTBEGIN(bugprone-use-after-move): checking the specified inactive state.
    CHECK_FALSE(contact.active());
    CHECK_FALSE(contact.geometry());
    // NOLINTEND(bugprone-use-after-move)
    CHECK(moved_contact.active());
}
TEST_CASE("Fixed nonplanar quad children and three incident bound faces retain shared sets without "
          "TE seams",
          "[UW-D04][UW-D07][UW-D08][UI-D03]") {
    core::ResourcePlan plan{budget};
    auto old = test_geometry::workspace(plan);
    auto work = workspace(plan);
    const std::array ids{s::PointId{1}, s::PointId{2}, s::PointId{3}, s::PointId{4}};
    const std::array bindings{test_geometry::binding(), test_geometry::binding(),
                              test_geometry::binding(), test_geometry::binding()};
    const std::array<double, 12> positions{0, 0, 0, 2, 0, 0, 2, 2, 1, 0, 2, 0};
    const std::array edges{s::EdgeInput{s::EdgeId{10}, {ids[0], ids[1]}},
                           s::EdgeInput{s::EdgeId{11}, {ids[1], ids[2]}},
                           s::EdgeInput{s::EdgeId{12}, {ids[2], ids[3]}},
                           s::EdgeInput{s::EdgeId{13}, {ids[3], ids[0]}},
                           s::EdgeInput{s::EdgeId{14}, {ids[0], ids[2]}}};
    const std::array faces{s::FaceInput{s::PhysicalFaceId{50}, s::BandId{100}}};
    const std::array bands{s::BandId{100}};
    const std::array quads{s::QuadInput{faces[0].id,
                                        ids,
                                        {edges[0].id, edges[1].id, edges[2].id, edges[3].id},
                                        edges[4].id,
                                        {s::FacetId{70}, s::FacetId{71}}}};
    const s::GeometryInput input{ids,
                                 bindings,
                                 s::SinglePositions{positions},
                                 s::SinglePositions{positions},
                                 edges,
                                 faces,
                                 bands,
                                 {},
                                 quads,
                                 s::SingleWord{0.25},
                                 test_geometry::context()};
    auto geometry =
        test_geometry::take(s::prepare_owned_geometry(input, test_geometry::limits(), *old, plan));
    const std::array roles{s::FaceRoleInput{faces[0].id, s::AeroFaceRole::Bound}};
    std::vector<s::EdgeBindingInput> edge_bindings;
    edge_bindings.reserve(edges.size());
    for (const auto& edge : edges)
        edge_bindings.push_back({edge.id, s::CanonicalEdge{}});
    auto topology =
        take(s::build_global_topology(geometry, {roles, edge_bindings, {}}, limits(), *work, plan));
    CHECK(take(s::check_global_geometry(topology, geometry, limits(), *work, plan)).active());
    const std::array diagonal{s::TrailingEdgeInput{edges[4].id, {ids[0], ids[2]}}};
    CHECK_FALSE(s::build_global_topology(geometry, {roles, edge_bindings, diagonal}, limits(),
                                         *work, plan));
    const std::array external{s::TrailingEdgeInput{edges[2].id, {ids[2], ids[3]}}};
    auto external_topology = take(s::build_global_topology(
        geometry, {roles, edge_bindings, external}, limits(), *work, plan));
    CHECK(external_topology.trailing_edges()[0].coefficient == -1);
    test_global_geometry::PairFixture triple{reference_case("edge")};
    triple.points.emplace_back(7);
    triple.bindings.push_back(test_geometry::binding());
    for (const auto value : {0.0, -2.0, -1.0}) {
        triple.high.push_back(value);
        triple.low.push_back(0);
    }
    triple.edges.emplace_back(s::EdgeId{30}, std::array{s::PointId{2}, s::PointId{7}});
    triple.edges.emplace_back(s::EdgeId{31}, std::array{s::PointId{7}, s::PointId{1}});
    triple.edge_bindings.push_back({s::EdgeId{30}, s::CanonicalEdge{}});
    triple.edge_bindings.push_back({s::EdgeId{31}, s::CanonicalEdge{}});
    std::vector<s::FaceInput> three_faces{triple.faces.begin(), triple.faces.end()};
    three_faces.emplace_back(s::PhysicalFaceId{52}, s::BandId{102});
    std::vector<s::BandId> three_bands{triple.bands.begin(), triple.bands.end()};
    three_bands.emplace_back(102);
    triple.triangles.emplace_back(s::FacetId{72}, s::PhysicalFaceId{52},
                                  std::array{s::PointId{1}, s::PointId{2}, s::PointId{7}},
                                  std::array{s::EdgeId{10}, s::EdgeId{30}, s::EdgeId{31}});
    const s::GeometryInput ti{triple.points,
                              triple.bindings,
                              s::DoublePositions{triple.high, triple.low},
                              s::DoublePositions{triple.high, triple.low},
                              triple.edges,
                              three_faces,
                              three_bands,
                              triple.triangles,
                              {},
                              s::SingleWord{0.25},
                              test_geometry::context()};
    auto shared =
        test_geometry::take(s::prepare_owned_geometry(ti, test_geometry::limits(), *old, plan));
    std::vector<s::FaceRoleInput> three_roles{triple.roles.begin(), triple.roles.end()};
    three_roles.push_back({s::PhysicalFaceId{52}, s::AeroFaceRole::Bound});
    auto shared_topology = take(s::build_global_topology(
        shared, {three_roles, triple.edge_bindings, {}}, limits(), *work, plan));
    CHECK(take(s::check_global_geometry(shared_topology, shared, limits(), *work, plan)).active());
    contact_budgets(shared_topology, shared, *work, 3, 3);
    const std::array seam{s::TrailingEdgeInput{s::EdgeId{10}, {s::PointId{1}, s::PointId{2}}}};
    CHECK_FALSE(s::build_global_topology(shared, {three_roles, triple.edge_bindings, seam},
                                         limits(), *work, plan));
}
TEST_CASE("Empty geometry has a qualified empty contact set and empty exact kinematics",
          "[UW-D07][N4-02][N4-03]") {
    core::ResourcePlan plan{budget};
    auto old = test_geometry::workspace(plan);
    auto work = workspace(plan);
    const s::GeometryInput input{
        {}, {}, s::SinglePositions{{}}, s::SinglePositions{{}},  {}, {}, {},
        {}, {}, s::SingleWord{0},       test_geometry::context()};
    auto source =
        test_geometry::take(s::prepare_owned_geometry(input, test_geometry::limits(), *old, plan));
    auto topology = take(s::build_global_topology(source, {{}, {}, {}}, limits(), *work, plan));
    CHECK(topology.trailing_edges().empty());
    CHECK(topology.trailing_point_count() == 0);
    auto strict = limits();
    strict.metadata = numerics::linear::RequireProvenBound{0};
    auto contacts = take(s::check_global_geometry(topology, source, strict, *work, plan));
    CHECK(take(contacts.work_plan()).candidates == 0);
    CHECK(take(contacts.work_plan()).traversal_executed_steps == 0);
    auto kinematics = take(s::check_alias_kinematics(source, {s::DofLayoutGeneration{0}, {}, {}},
                                                     strict, *work, plan));
    CHECK(kinematics.active());
    CHECK_FALSE(kinematics.velocity(s::PointId{0}, 0));
}
TEST_CASE("Original word extremes preserve nondegeneracy and exact velocity beyond binary64 range",
          "[UW-D07][N4-02][N4-06]") {
    for (const auto magnitude :
         {std::numeric_limits<double>::max(), std::numeric_limits<double>::denorm_min()}) {
        INFO("source word " << magnitude);
        test_geometry::TriangleFixture fixture;
        fixture.two_word = true;
        fixture.reference = {0, 0, 0, magnitude, 0, 0, 0, magnitude, 0};
        fixture.current = fixture.reference;
        fixture.reference_low = fixture.reference;
        fixture.current_low = fixture.current;
        core::ResourcePlan plan{budget};
        auto old = test_geometry::workspace(plan);
        auto work = workspace(plan);
        auto source = test_geometry::take(
            s::prepare_owned_geometry(fixture.input(), test_geometry::limits(), *old, plan));
        const std::array roles{s::FaceRoleInput{fixture.faces[0].id, s::AeroFaceRole::Bound}};
        const std::array edges{s::EdgeBindingInput{fixture.edges[0].id, s::CanonicalEdge{}},
                               s::EdgeBindingInput{fixture.edges[1].id, s::CanonicalEdge{}},
                               s::EdgeBindingInput{fixture.edges[2].id, s::CanonicalEdge{}}};
        auto topology =
            take(s::build_global_topology(source, {roles, edges, {}}, limits(), *work, plan));
        CHECK(take(s::check_global_geometry(topology, source, limits(), *work, plan)).active());
        const std::array<s::ScalarWords, 3> axis{s::SingleWord{1}, s::SingleWord{0},
                                                 s::SingleWord{0}};
        const std::array dofs{s::GlobalDofInput{core::DofId{1},
                                                {1, 0, 0, 0, 0, 0, 0},
                                                s::DofBasisId{1},
                                                axis,
                                                s::DoubleWord{magnitude, magnitude}}};
        const std::array entries{
            s::GlobalJacobianEntry{core::DofId{1}, s::DoubleWord{magnitude, magnitude}}};
        const std::array<s::ScalarWords, 3> zero{s::SingleWord{0}, s::SingleWord{0},
                                                 s::SingleWord{0}};
        const std::array points{
            s::PointKinematicsInput{fixture.points[0], {entries, {}, {}}, zero},
            s::PointKinematicsInput{fixture.points[1], {entries, {}, {}}, zero},
            s::PointKinematicsInput{fixture.points[2], {entries, {}, {}}, zero}};
        auto result = take(s::check_alias_kinematics(
            source, {s::DofLayoutGeneration{0}, dofs, points}, limits(), *work, plan));
        const auto velocity = take(result.velocity(fixture.points[0], 0));
        boost::multiprecision::cpp_int actual = 0;
        for (auto i = velocity.limbs.size(); i > 0; --i) {
            actual <<= 32U;
            actual += velocity.limbs[i - 1];
        }
        actual <<= 32U * velocity.first_limb;
        // max = (2^53-1)*2^971; each supplied two-word sum is twice max.
        const boost::multiprecision::cpp_int expected =
            magnitude == std::numeric_limits<double>::max()
                ? boost::multiprecision::cpp_int{((boost::multiprecision::cpp_int{1} << 53U) - 1) *
                                                 ((boost::multiprecision::cpp_int{1} << 53U) - 1)}
                      << 6240U
                : boost::multiprecision::cpp_int{1} << 2150U;
        CHECK(velocity.sign == 1);
        CHECK(actual == expected);
    }
}
