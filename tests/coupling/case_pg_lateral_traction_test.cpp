#include "../io/reference_declarations_fixture.hpp"
#include "../structure/pg_partition_test_support.hpp"
#include "coupling/case_pg_lateral_traction.hpp"
#include "structure/exact_polynomial_detail.hpp"

#include <boost/multiprecision/cpp_int.hpp>
#include <catch2/catch_test_macros.hpp>

#include <limits>
#include <utility>

namespace {
using namespace uvlas;
using namespace uvlas::coupling::detail;
namespace fixture = test::declaration_fixture;
namespace pg = structure;
using Q = boost::multiprecision::cpp_rational;
using Vector = std::array<Q, 3>;
using Code = core::FailureCode;

const Q thickness{127, 50'000};
const Q axial_traction{102550048828125LL, 524288};
const Q tolerance{1, boost::multiprecision::cpp_int{"1000000000000000000"}};

Q exact(const pg::ExactRational& value) {
    const auto& backend = pg::detail::ExactFactory::value(value).backend();
    return {boost::multiprecision::cpp_int{backend.num()},
            boost::multiprecision::cpp_int{backend.denom()}};
}

struct Arithmetic {
    pg::PgReferenceLimits limits{partition_test::exact_limits(), 100'000'000, 20'000'000, {}};
    core::ResourcePlan memory{fixture::budget};
    core::Workspace workspace{core::PlannedBuffer<std::byte>{}};

    Arithmetic() {
        limits.algebra.maximum_value_bits = 16'384;
        limits.algebra.maximum_work_bits = 65'536;
    }

    pg::ExactRational own(const Q& q) {
        return fixture::take(pg::ExactRational::from_integer_pair(
            numerator(q).str(), denominator(q).str(), limits.algebra, memory));
    }

    io::SiCase si(const fixture::JsonToken& input) {
        const auto checked = fixture::checked(input);
        std::array<std::uint64_t, static_cast<std::size_t>(io::SiQuota::Count)> quotas{};
        quotas.fill(fixture::budget);
        const io::SiLimits caps{quotas, limits.algebra};
        const auto catalog = fixture::take(io::make_si_catalog(checked, caps, memory, nullptr));
        std::vector<io::SiRequirementRow> requirements;
        for (std::uint64_t i = 0; i < fixture::take(catalog.size()); ++i) {
            const auto entry = fixture::take(catalog.entry(i));
            requirements.push_back({entry.key,
                                    {io::SiDecimalInput{"1"}, entry.recipe.dimensions},
                                    io::SiDecimalInput{"1e-25"},
                                    io::SiDecimalInput{"1e-25"}});
        }
        const auto request =
            fixture::take(io::make_si_request(catalog, requirements, caps, memory, nullptr));
        return fixture::take(
            io::convert_case_to_si(checked, catalog, request, caps, memory, nullptr));
    }

    CasePgChartLimits chart_limits() const {
        constexpr auto profile = std::uint64_t{48} << 30U;
        constexpr std::uint64_t steps = 10'000'000'000'000;
        const CasePgReferenceLimits geometry{{limits.algebra, 4, {}},
                                             {limits.algebra, 2000, 20'000'000, 100'000, {}},
                                             partition_test::limits(),
                                             limits,
                                             100'000,
                                             {}};
        auto factor = limits;
        factor.algebra.maximum_value_bits = 65'536;
        factor.algebra.maximum_work_bits = 262'144;
        return {geometry,
                {profile, steps, 128, {}},
                {768, 704},
                {factor, steps},
                {768, 192},
                {profile, steps, {}},
                factor,
                {limits.algebra, 10'000, case_motion_frame_interval_scratch_bytes},
                {factor, 120'000},
                {factor, 1, 1, 500, 120'000},
                {factor, 512, 960, 140'000},
                100'000};
    }

    core::Result<CaseComponentPgChart> chart(const io::SiCase& input, bool different_rule = false) {
        const auto sixth = own(Q{1, 6}), third = own(Q{1, 3});
        std::array<pg::ExactRational, 6> positions;
        const std::array<pg::PgLineQuadratureNode, 6> plane{{{positions[0], sixth},
                                                             {positions[1], sixth},
                                                             {positions[2], sixth},
                                                             {positions[3], sixth},
                                                             {positions[4], sixth},
                                                             {positions[5], sixth}}};
        for (std::size_t i = 0; i < positions.size(); ++i)
            positions[i] = own(different_rule ? Q{2 * i + 1, 13} : Q{2 * i + 1, 12});
        const auto quarter = own(Q{1, 4}), half = own(Q{1, 2}), upper = own(Q{3, 4});
        const std::array depth{pg::PgLineQuadratureNode{quarter, third},
                               pg::PgLineQuadratureNode{half, third},
                               pg::PgLineQuadratureNode{upper, third}};
        const std::array rules{CasePgRuleBinding{1, pg::PgVolumeProductRule{plane, plane, depth}}};
        return prepare_case_component_pg_chart(input, rules, chart_limits(), workspace, memory);
    }
};

fixture::JsonToken rectangle() {
    auto input = fixture::root("structural_static");
    auto& model = input.field("model");
    auto& geometry = model.field("geometry").elements.front();
    geometry.field("frame").set("origin", fixture::json("[1,2,3]"));
    geometry.field("frame").set("rotation", fixture::json("[[0,0,1],[1,0,0],[0,1,0]]"));
    geometry.set("reference_length", fixture::json("2"));
    geometry.field("thickness").set("values", fixture::json("[0.00254,0.00254,0.00254,0.00254]"));
    geometry.field("boundaries").elements[1].set("vertices", fixture::json(R"([
        [{"numerator":"1","denominator":"1"},{"numerator":"1","denominator":"1"}],
        [{"numerator":"1","denominator":"1"},{"numerator":"1","denominator":"5"}],
        [{"numerator":"1","denominator":"1"},{"numerator":"0","denominator":"1"}]
    ])"));
    auto& layers = model.field("shells").elements.front().field("layers").elements;
    auto upper = layers.front();
    layers.front().set("id", fixture::string("11"));
    layers.front().set("interval", fixture::json(R"([
        {"numerator":"-1","denominator":"1"},{"numerator":"-1","denominator":"3"}
    ])"));
    upper.set("id", fixture::string("29"));
    upper.set("interval", fixture::json(R"([
        {"numerator":"-1","denominator":"3"},{"numerator":"1","denominator":"1"}
    ])"));
    layers.push_back(std::move(upper));

    auto load = fixture::fragment("Load.traction");
    load.set("id", fixture::string("31"));
    auto& region = load.field("region");
    region.set("kind", fixture::string("edge"));
    region.set("side", fixture::string("interior"));
    region.set("boundary", fixture::json("1"));
    load.field("spatial").set("value", fixture::json("[0,195598695.4271793365478515625,0]"));
    auto second = load;
    second.set("id", fixture::string("17"));
    second.field("spatial").set("value", fixture::json("[8,3,-5]"));
    second.field("multiplier").set("value", fixture::json("[-2]"));
    model.field("shells").elements.front().field("loads").elements.push_back(std::move(load));
    input.field("analysis").field("loads").elements.push_back(std::move(second));
    return input;
}

Q power(Q value, core::Index n) {
    Q product = 1;
    for (core::Index i = 0; i < n; ++i)
        product *= value;
    return product;
}

// Independently integrate the cardinal polynomials of the adopted 22 rational nodes.
struct EdgeRule {
    static constexpr std::array numerators{
        0,      5856,   23293,  51921,  91101,  139958, 197400, 262144, 332744,  407623,  485108,
        563468, 640953, 715832, 786432, 851176, 908618, 957475, 996655, 1025283, 1042720, 1048576};
    std::array<Q, numerators.size()> u{}, w{};
    std::array<pg::ExactRational, numerators.size()> nodes{}, weights{};
    std::vector<pg::PgLineQuadratureNode> rule;

    explicit EdgeRule(Arithmetic& a) {
        for (std::size_t i = 0; i < u.size(); ++i)
            u[i] = Q{numerators[i], 1 << 20};
        for (std::size_t i = 0; i < u.size(); ++i) {
            std::array<Q, numerators.size()> cardinal{};
            cardinal[0] = 1;
            std::size_t degree = 0;
            for (std::size_t j = 0; j < u.size(); ++j) {
                if (i == j)
                    continue;
                const Q divisor = u[i] - u[j];
                for (std::size_t n = ++degree; n > 0; --n)
                    cardinal[n] = (cardinal[n - 1] - u[j] * cardinal[n]) / divisor;
                cardinal[0] *= -u[j] / divisor;
            }
            for (std::size_t n = 0; n <= degree; ++n)
                w[i] += cardinal[n] / (n + 1);
            REQUIRE(w[i] > 0);
            nodes[i] = a.own(u[i]);
            weights[i] = a.own(w[i]);
            rule.push_back({nodes[i], weights[i]});
        }
        for (core::Index degree = 0; std::cmp_less(degree, u.size()); ++degree) {
            Q moment;
            for (std::size_t i = 0; i < u.size(); ++i)
                moment += w[i] * power(u[i], degree);
            CHECK(moment == Q{1, degree + 1});
        }
    }
};

struct Interval {
    Q lower{}, upper{};

    void add(const pg::PgMaterialEnclosure& value, const Q& direction) {
        lower += direction * exact(direction < 0 ? value.upper : value.lower);
        upper += direction * exact(direction < 0 ? value.lower : value.upper);
    }

    void check(const Q& expected) const {
        INFO("independent: " << expected << " interval: [" << lower << ", " << upper << ']');
        CHECK(lower <= expected);
        CHECK(expected <= upper);
        CHECK(upper - lower <= tolerance);
    }
};

Vector cross(const Vector& a, const Vector& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

struct IndependentColumns {
    std::vector<Q> unit_pascals, depth_direction;
};

IndependentColumns integrate_polynomials(const pg::PgPreparedPatch& patch, Arithmetic& a) {
    const auto* basis = fixture::take(patch.conditioning().reference_basis());
    const auto* rank = fixture::take(basis->rank());
    const auto rows = fixture::take(rank->rows());
    const auto raw = fixture::take(rank->provenance());
    std::vector<Q> integrals(patch.translation().size()), depth(integrals.size());
    const Q z_scale = thickness / 4;
    for (std::size_t i = 0; i < integrals.size(); ++i) {
        const pg::PgRankBasisColumn column{static_cast<core::Index>(i)};
        for (const auto k : fixture::take(rank->basis_indices(column))) {
            const auto& row = rows[k];
            if (row.cell.first != 1)
                continue;
            REQUIRE(row.powers.second <= 21);
            if (row.powers.second % 2 != 0 || row.field == pg::PgCoefficientField::P)
                continue;
            // On this flat rectangle, u=1, dt=dv/4 and dA=H dt dζ. The
            // two t-cells are separate [-1,1] intervals in their own v.
            const Q thickness_integral =
                row.field == pg::PgCoefficientField::C ? Q{z_scale * z_scale / 3} : Q{2};
            const Q plane_integral{1, 2 * (row.powers.second + 1)};
            const auto coefficient =
                fixture::take(rank->basis(column, k, a.limits.algebra, a.memory));
            integrals[i] += exact(coefficient) * thickness * plane_integral * thickness_integral;
        }
        // The tensor C hats form a partition of unity, so (2/z_scale^2) times
        // their sum in the unconditioned basis is the virtual displacement ζ².
        for (std::size_t j = 0; j < raw.size(); ++j)
            if (raw[j].kind == pg::PgCoefficientKind::BasicC)
                depth[i] += exact(fixture::take(rank->original(
                    column, {static_cast<core::Index>(j)}, a.limits.algebra, a.memory)));
        depth[i] *= 2 / (z_scale * z_scale);
    }
    const auto transform = [&](std::size_t i, std::size_t j) {
        return exact(fixture::take(patch.conditioning().transform_coefficient(
            {static_cast<core::Index>(i)}, {static_cast<core::Index>(j)}, a.limits.algebra,
            a.memory)));
    };
    std::vector<Q> conditioned(integrals.size());
    for (std::size_t j = 0; j < integrals.size(); ++j)
        for (std::size_t i = 0; i <= j; ++i)
            conditioned[j] += transform(i, j) * integrals[i];
    for (std::size_t end = depth.size(); end > 0; --end) {
        const auto i = end - 1;
        for (std::size_t j = i + 1; j < depth.size(); ++j)
            depth[i] -= transform(i, j) * depth[j];
        depth[i] /= transform(i, i);
    }
    return {std::move(conditioned), std::move(depth)};
}

void check_provenance(std::span<const CasePgLateralTractionPoint> points, const EdgeRule& edge,
                      const pg::PgPreparedPatch& patch, Arithmetic& a) {
    const std::array range_start{Q{1}, Q{1, 2}, Q{1, 5}};
    const std::array range_end{Q{1, 2}, Q{1, 5}, Q{0}};
    const std::array depths{Q{0}, Q{1, 2}, Q{1}};
    const std::array depth_weights{Q{1, 6}, Q{2, 3}, Q{1, 6}};
    const std::array layer_lower{Q{-1}, Q{-1, 3}};
    const std::array layer_width{Q{2, 3}, Q{4, 3}};
    REQUIRE(points.size() == 2 * 3 * 2 * edge.rule.size() * depths.size());
    std::size_t index = 0;
    const auto basis = fixture::take(patch.conditioning().reference_basis());
    for (const auto id : {31U, 17U})
        for (std::size_t atom = 0; atom < 3; ++atom)
            for (std::size_t layer = 0; layer < 2; ++layer)
                for (std::size_t u = 0; u < edge.rule.size(); ++u)
                    for (std::size_t v = 0; v < depths.size(); ++v) {
                        const auto& point = points[index++];
                        CHECK(point.load == id);
                        CHECK(point.boundary == 1);
                        CHECK(point.segment == (atom == 2 ? 1 : 0));
                        CHECK(point.layer == core::LayerId{layer == 0 ? 11U : 29U});
                        CHECK(exact(point.source_along) == edge.u[u]);
                        CHECK(exact(point.face_along) == 1 - edge.u[u]);
                        CHECK(exact(point.depth) == depths[v]);
                        CHECK(exact(point.weight) == edge.w[u] * depth_weights[v]);
                        if (id != 31 || u != 6 || v != 1)
                            continue;
                        const auto location = fixture::take(basis->evaluate_geometry(
                            pg::PgReferenceLateralPoint{point.face, point.face_along, point.depth},
                            a.limits, a.memory));
                        const Q t =
                            range_start[atom] + edge.u[u] * (range_end[atom] - range_start[atom]);
                        const Q zeta = layer_lower[layer] + depths[v] * layer_width[layer];
                        CHECK(exact(location.position_metres[0].rational) +
                                  exact(location.position_metres[0].radical) ==
                              1 + thickness * zeta / 2);
                        CHECK(exact(location.position_metres[1].rational) +
                                  exact(location.position_metres[1].radical) ==
                              4);
                        CHECK(exact(location.position_metres[2].rational) +
                                  exact(location.position_metres[2].radical) ==
                              3 + 2 * t);
                        CHECK(exact(location.g) == 1);
                    }
}

void check_action(const CasePgLateralTractionForce& force, const pg::PgPreparedPatch& patch,
                  Arithmetic& a) {
    const auto polynomial = integrate_polynomials(patch, a);
    const auto rows = force.newtons.values();
    REQUIRE(rows.size() == polynomial.unit_pascals.size());
    const Vector traction{-16, axial_traction - 6, 10};
    const Vector origin{2, 1, 3}, shift{Q{1, 2}, Q{-3, 4}, 1}, rotation{Q{1, 4}, Q{-1, 2}, Q{3, 4}};
    std::array<Interval, 3> sum{}, moment{};
    Interval rigid, affine, square, extension;
    for (std::size_t i = 0; i < rows.size(); ++i) {
        CAPTURE(i);
        const Q translation = exact(patch.translation()[i]);
        Vector radius{};
        for (std::size_t axis = 0; axis < 3; ++axis)
            radius[axis] =
                exact(patch.reference_coordinates()[i][axis]) - origin[axis] * translation;
        const auto displacement = cross(rotation, radius);
        for (std::size_t axis = 0; axis < 3; ++axis) {
            CAPTURE(axis);
            Interval original;
            original.add(rows[i][axis], 1);
            original.check(polynomial.unit_pascals[i] * traction[axis]);
            sum[axis].add(rows[i][axis], translation);
            const auto j = (axis + 1) % 3, k = (axis + 2) % 3;
            moment[axis].add(rows[i][k], radius[j]);
            moment[axis].add(rows[i][j], -radius[k]);
            rigid.add(rows[i][axis], shift[axis] * translation + displacement[axis]);
        }
        affine.add(rows[i][1], (exact(patch.reference_coordinates()[i][1]) - 2 * translation) / 2);
        square.add(rows[i][1], polynomial.depth_direction[i]);
        extension.add(rows[i][1], translation / 512);
    }
    const Q area = 2 * thickness;
    const Vector expected_force{area * traction[0], area * traction[1], area * traction[2]};
    const Vector expected_moment = cross(Vector{-1, 3, 1}, expected_force);
    Q expected_rigid;
    for (std::size_t i = 0; i < 3; ++i) {
        sum[i].check(expected_force[i]);
        moment[i].check(expected_moment[i]);
        expected_rigid += shift[i] * expected_force[i] + rotation[i] * expected_moment[i];
    }
    CHECK(area * axial_traction == Q{4167633984375LL, 4194304});
    rigid.check(expected_rigid);
    affine.check(expected_force[1]);
    square.check(expected_force[1] / 3);
    extension.check(expected_force[1] / 512);
}

TEST_CASE("Case edge tractions integrate the original PG columns over all side layers",
          "[FD-05][BC-01][LOAD-01][V-LOAD][integration][verification][negative]") {
    Arithmetic a;
    const auto input = rectangle();
    const auto si = a.si(input);
    const auto owner = fixture::take(a.chart(si));
    const auto bound = fixture::take(owner.point_patch(0));
    const auto zero = a.own(0), half = a.own(Q{1, 2}), one = a.own(1);
    const auto sixth = a.own(Q{1, 6}), two_thirds = a.own(Q{2, 3});
    const std::array depth{pg::PgLineQuadratureNode{zero, sixth},
                           pg::PgLineQuadratureNode{half, two_thirds},
                           pg::PgLineQuadratureNode{one, sixth}};
    const EdgeRule edge{a};
    const CasePgLateralTractionRule rule{edge.rule, depth};
    const CasePgLateralTractionLimits limits{a.limits, 2, 6, 2, 25, 792};
    CasePgLateralTractionDiagnostic diagnostic;
    auto prepared = CasePgLateralTraction::prepare(si, bound, rule, limits, a.memory, &diagnostic);
    REQUIRE(prepared);
    CHECK_FALSE(diagnostic.load);
    check_provenance(prepared->points(), edge, bound.patch(), a);
    const auto width = a.own(tolerance);
    const CasePgLateralTractionAccuracy accuracy{160, width};
    const auto force = prepared->evaluate(bound, 4.25, accuracy, limits, a.memory, &diagnostic);
    REQUIRE(force);
    CHECK(force->patch == fixture::take(bound.patch().info()).identity);
    CHECK(force->shell == numerics::PatchId{1});
    CHECK_FALSE(diagnostic.load);
    check_action(*force, bound.patch(), a);

    const auto exhausted = [&](const CasePgLateralTractionLimits& capped) {
        const auto failed = CasePgLateralTraction::prepare(si, bound, rule, capped, a.memory);
        REQUIRE_FALSE(failed);
        CHECK(failed.error().code == Code::ResourceExceeded);
    };
    auto capped = limits;
    --capped.maximum_loads;
    exhausted(capped);
    capped = limits;
    --capped.maximum_atoms;
    exhausted(capped);
    capped = limits;
    --capped.maximum_layers;
    exhausted(capped);
    capped = limits;
    --capped.maximum_rule_nodes;
    exhausted(capped);
    capped = limits;
    --capped.maximum_points;
    exhausted(capped);
    capped = limits;
    capped.reference.maximum_scalar_operations = 0;
    exhausted(capped);
    capped = limits;
    capped.reference.maximum_partition_visits = 0;
    exhausted(capped);

    const auto check_failure = [&](auto&& failed, Code expected) {
        REQUIRE_FALSE(failed);
        CHECK(failed.error().code == expected);
    };
    core::ResourcePlan tiny{0};
    check_failure(CasePgLateralTraction::prepare(si, bound, rule, limits, tiny),
                  Code::ResourceExceeded);
    check_failure(prepared->evaluate(bound, 0, accuracy, limits, tiny), Code::ResourceExceeded);
    const auto stopped = partition_test::requested_stop_token();
    auto cancelled = limits;
    cancelled.reference.cancellation = stopped;
    check_failure(CasePgLateralTraction::prepare(si, bound, rule, cancelled, a.memory),
                  Code::Cancelled);
    check_failure(prepared->evaluate(bound, 0, accuracy, cancelled, a.memory), Code::Cancelled);
    check_failure(prepared->evaluate(bound, std::numeric_limits<double>::infinity(), accuracy,
                                     limits, a.memory),
                  Code::Nonfinite);
    check_failure(prepared->evaluate(bound, std::numeric_limits<double>::quiet_NaN(), accuracy,
                                     limits, a.memory),
                  Code::Nonfinite);
    check_failure(prepared->evaluate(bound, 0, {32768, width}, limits, a.memory),
                  Code::PrecisionUnsupported);

    for (const std::array invalid :
         {std::array{pg::PgLineQuadratureNode{half, one}, pg::PgLineQuadratureNode{zero, one}},
          std::array{pg::PgLineQuadratureNode{half, half}, pg::PgLineQuadratureNode{half, half}},
          std::array{pg::PgLineQuadratureNode{zero, zero}, pg::PgLineQuadratureNode{one, one}},
          std::array{pg::PgLineQuadratureNode{zero, half}, pg::PgLineQuadratureNode{one, sixth}}})
        check_failure(CasePgLateralTraction::prepare(si, bound, {invalid, depth}, limits, a.memory),
                      Code::InvalidInput);
    check_failure(CasePgLateralTraction::prepare(si, bound, {{}, depth}, limits, a.memory),
                  Code::InvalidInput);

    const auto changed_basis = fixture::take(a.chart(si, true));
    const auto changed_bound = fixture::take(changed_basis.point_patch(0));
    CHECK(changed_bound.patch().translation().size() == bound.patch().translation().size());
    check_failure(prepared->evaluate(changed_bound, 0, accuracy, limits, a.memory),
                  Code::StaleGeneration);

    auto foreign_input = input;
    foreign_input.field("model")
        .field("shells")
        .elements.front()
        .field("loads")
        .elements.front()
        .set("area", fixture::string("current"));
    const auto foreign_si = a.si(foreign_input);
    check_failure(CasePgLateralTraction::prepare(foreign_si, bound, rule, limits, a.memory),
                  Code::StaleGeneration);
    const auto foreign_owner = fixture::take(a.chart(foreign_si));
    const auto foreign_bound = fixture::take(foreign_owner.point_patch(0));
    check_failure(prepared->evaluate(foreign_bound, 0, accuracy, limits, a.memory),
                  Code::StaleGeneration);
    check_failure(CasePgLateralTraction::prepare(foreign_si, foreign_bound, rule, limits, a.memory,
                                                 &diagnostic),
                  Code::Unsupported);
    CHECK(diagnostic.load == 31);
    CHECK(diagnostic.boundary == 1);

    // Owned partitions reject a named line with incomplete physical coverage
    // or positive self-overlap before a subsequent traction can bind to it.
    for (std::string_view vertices : {
             R"([[{"numerator":"1","denominator":"1"},{"numerator":"1","denominator":"1"}],
                 [{"numerator":"3","denominator":"4"},{"numerator":"1","denominator":"5"}]])",
             R"([[{"numerator":"1","denominator":"1"},{"numerator":"1","denominator":"1"}],
                 [{"numerator":"1","denominator":"1"},{"numerator":"1","denominator":"5"}],
                 [{"numerator":"1","denominator":"1"},{"numerator":"1","denominator":"2"}]])"}) {
        auto bad_case = input;
        bad_case.field("model")
            .field("geometry")
            .elements.front()
            .field("boundaries")
            .elements[1]
            .set("vertices", fixture::json(vertices));
        check_failure(a.chart(a.si(bad_case)), Code::GeometryInvalid);
    }

    const auto third = a.own(Q{1, 3});
    const std::array asymmetric{pg::PgLineQuadratureNode{third, one}};
    const auto narrow = fixture::take(
        CasePgLateralTraction::prepare(si, bound, {asymmetric, asymmetric}, limits, a.memory));
    REQUIRE(narrow.points().size() == 12);
    CHECK(exact(narrow.points().front().face_along) == Q{2, 3});
    check_failure(narrow.evaluate(bound, 0, {0, zero}, limits, a.memory),
                  Code::PrecisionUnresolved);
    capped = limits;
    capped.reference.maximum_partition_visits = force->newtons.size();
    check_failure(narrow.evaluate(bound, 0, accuracy, capped, a.memory, &diagnostic),
                  Code::ResourceExceeded);
    CHECK(diagnostic.load == 31);
    CHECK(diagnostic.layer == core::LayerId{11});
    REQUIRE(narrow.evaluate(bound, -0.75, accuracy, limits, a.memory));

    // The original accepted preparation remains readable after each rejected attempt.
    CHECK(prepared->points().size() == 792);
    const auto saved = std::move(*prepared);
    check_failure(prepared->evaluate(bound, 0, accuracy, limits, a.memory), Code::StaleGeneration);
    CHECK(saved.points().size() == 792);
}

TEST_CASE("Case edge tractions enclose a curved reference finite rule with an irrational area",
          "[FD-05][V-LOAD][integration][verification]") {
    Arithmetic a;
    auto input = rectangle();
    auto& geometry = input.field("model").field("geometry").elements.front();
    geometry.set("degree", fixture::json("[0,2]"));
    geometry.set("coefficients", fixture::json(R"([
        {"powers":[0,2],"coefficient":{"numerator":"1","denominator":"2"}}
    ])"));
    input.field("analysis").field("loads").elements.clear();
    const auto si = a.si(input);
    const auto owner = fixture::take(a.chart(si));
    const auto bound = fixture::take(owner.point_patch(0));
    const auto third = a.own(Q{1, 3}), one = a.own(1), width = a.own(tolerance);
    const std::array nodes{pg::PgLineQuadratureNode{third, one}};
    const CasePgLateralTractionLimits limits{a.limits, 1, 3, 2, 2, 6};
    const auto prepared =
        fixture::take(CasePgLateralTraction::prepare(si, bound, {nodes, nodes}, limits, a.memory));
    REQUIRE(prepared.points().size() == 6);
    const auto force = fixture::take(prepared.evaluate(bound, 0, {160, width}, limits, a.memory));

    // For X = 2(s,t,t²/2) + Hζ n/2 before the rigid frame,
    // n = (0,-t,1)/sqrt(g) with g = 1+t², and n_t = -(0,1,t)/g^(3/2).
    // These come from Cartesian differentiation, independently of PG point evaluation.
    const std::array edge_length{Q{1, 2}, Q{3, 10}, Q{1, 5}};
    const std::array t{Q{5, 6}, Q{2, 5}, Q{2, 15}};
    const std::array layers{Q{2, 3}, Q{4, 3}};
    const std::array zeta{Q{-7, 9}, Q{1, 9}};
    const boost::multiprecision::cpp_int scale = boost::multiprecision::cpp_int{1} << 192;
    Interval expected;
    for (std::size_t i = 0; i < t.size(); ++i) {
        const Q g = 1 + t[i] * t[i];
        const boost::multiprecision::cpp_int square = numerator(g) * scale * scale / denominator(g);
        const boost::multiprecision::cpp_int low = sqrt(square);
        for (std::size_t j = 0; j < layers.size(); ++j) {
            const Q factor = axial_traction * thickness * edge_length[i] * layers[j] / 2;
            const Q curvature = thickness * zeta[j] / (2 * g);
            expected.lower += factor * (2 * Q{low, scale} - curvature);
            expected.upper += factor * (2 * Q{low + 1, scale} - curvature);
        }
    }
    Interval measured;
    for (std::size_t i = 0; i < force.newtons.size(); ++i)
        measured.add(force.newtons.values()[i][1], exact(bound.patch().translation()[i]));
    CHECK(measured.lower <= expected.lower);
    CHECK(expected.upper <= measured.upper);
    CHECK(measured.upper - measured.lower <= tolerance);
}

} // namespace
