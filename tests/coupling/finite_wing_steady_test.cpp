#include "aero/geometry_snapshot_fixture.hpp"
#include "coupling/case_reference_graph.hpp"
#include "coupling/case_uniform_flow.hpp"
#include "coupling/case_uniform_flow_projection.hpp"
#include "coupling/discrete_nonpenetration.hpp"
#include "io/case_schema_test_adapter.hpp"
#include "io/si_rational_projection.hpp"

#include <boost/multiprecision/cpp_int.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <numbers>
#include <ranges>
#include <tuple>
#include <uvlas/numerics/dense_solve.hpp>
#include <vector>

namespace {
using namespace uvlas;
using namespace uvlas::test::schema_adapter;
namespace a = aero::detail;
namespace p = a::pchip;
namespace cf = coupling::detail;

constexpr std::uint64_t budget = std::uint64_t{1} << 48U;
const p::PchipInputStamp stamp{p::PchipIssuerId{71}, p::PchipInputGeneration{8}};

structure::ExactLimits exact_limits() noexcept {
    return {100000,
            100000,
            100000,
            64,
            64,
            100000,
            200000,
            structure::ExactEstimateAdmission{budget},
            structure::ExactTimingMode::Disabled};
}
const cf::CaseUniformFlowLimits limits{
    exact_limits(),
    {{128, 100000, budget, 10000}, numerics::linear::EstimateAdmission{budget}},
    100};

template <class T>
T take(core::Result<T> value) {
    if (!value)
        FAIL("Unexpected case-flow failure " << static_cast<int>(value.error().code));
    return std::move(*value);
}

template <class T, class Diagnostic>
T take_io(core::Result<T> value, const Diagnostic& diag) {
    if (!value) {
        const auto path = diag.path();
        if (path)
            FAIL("Unexpected case-flow input at " << *path << " code "
                                                  << static_cast<int>(value.error().code));
        FAIL("Unexpected case-flow preflight " << static_cast<int>(value.error().code));
    }
    return std::move(*value);
}

io::SiCase materialize(const JsonToken& source, io::SourceId id) {
    core::ResourcePlan memory{budget};
    const auto wire = source.wire();
    const io::InputLimits syntax{
        1000000, 256,     100000,  100000,
        1000000, 1000000, 1000000, {io::ParserAllocationMode::EstimateAdmission, budget}};
    std::optional<io::InputDiagnostic> parsing;
    const auto doc = take(io::read_input_document(
        std::as_bytes(std::span{wire.data(), wire.size()}), id, syntax, memory, parsing));
    std::array<std::uint64_t, static_cast<std::size_t>(io::SchemaQuota::Count)> schema_caps{};
    schema_caps.fill(budget);
    io::SchemaDiagnostic schema_diagnostic;
    const auto checked = take_io(
        io::validate_case_schema(doc, io::SchemaLimits{schema_caps}, memory, &schema_diagnostic),
        schema_diagnostic);

    std::array<std::uint64_t, static_cast<std::size_t>(io::SiQuota::Count)> si_caps{};
    si_caps.fill(budget);
    const io::SiLimits si_limits{si_caps, exact_limits()};
    io::SiDiagnostic si_diagnostic;
    const auto catalog =
        take_io(io::make_si_catalog(checked, si_limits, memory, &si_diagnostic), si_diagnostic);
    std::vector<io::SiRequirementRow> requirements;
    for (std::uint64_t i = 0; i < take(catalog.size()); ++i) {
        const auto quantity = take(catalog.entry(i));
        requirements.push_back({quantity.key,
                                {io::SiDecimalInput{"1"}, quantity.recipe.dimensions},
                                io::SiDecimalInput{"1e-25"},
                                io::SiDecimalInput{"1e-25"}});
    }
    const auto request =
        take_io(io::make_si_request(catalog, requirements, si_limits, memory, &si_diagnostic),
                si_diagnostic);
    return take_io(
        io::convert_case_to_si(checked, catalog, request, si_limits, memory, &si_diagnostic),
        si_diagnostic);
}

struct Number {
    Number(std::uint32_t numerator, std::uint32_t denominator = 1, int sign = 1) noexcept
        : n{numerator}, d{denominator}, sign_(numerator == 0 ? 0 : sign) {}
    [[nodiscard]] p::RationalView view() const& noexcept {
        return {sign_, sign_ == 0 ? std::span<const std::uint32_t>{} : n, d};
    }
    [[nodiscard]] p::RationalView view() const&& = delete;
    std::array<std::uint32_t, 1> n, d;
    int sign_;
};

using Big = boost::multiprecision::cpp_int;
using Rat = boost::multiprecision::cpp_rational;
Rat two_power(std::int64_t exponent) {
    if (exponent >= 0)
        return Rat{Big{1} << static_cast<unsigned>(exponent)};
    return Rat{1} / (Big{1} << static_cast<unsigned>(-exponent));
}

/// Independent exact reconstruction of a binary word using its binary decomposition.
Rat rational_word(double word) {
    REQUIRE(std::isfinite(word));
    int exponent = 0;
    const auto fraction = std::frexp(word, &exponent);
    const auto significand = static_cast<std::int64_t>(std::ldexp(fraction, 53));
    return Rat{significand} * two_power(exponent - 53);
}

template <class T>
T take_global(a::GlobalResult<T> result) {
    if (!result)
        FAIL("Unexpected global qualification failure "
             << static_cast<int>(result.error().origin.cause.code));
    return std::move(*result);
}

struct FiniteWingMesh {
    explicit FiniteWingMesh(const std::array<std::array<double, 3>, 20>& coordinates) {
        const auto point = [](std::size_t chord, std::size_t span) {
            return a::PointId{1 + 5 * chord + span};
        };
        for (std::size_t i = 0; i < coordinates.size(); ++i) {
            ids.emplace_back(i + 1);
            bindings.push_back(test_geometry::binding());
            const auto x = coordinates[i];
            xyz.insert(xyz.end(), x.begin(), x.end());
        }
        std::map<std::array<a::PointId, 2>, a::EdgeId> by_ends;
        const auto edge = [&](a::PointId from, a::PointId to) {
            const std::array key{std::min(from, to), std::max(from, to)};
            const a::EdgeId new_id{1000 + edges.size()};
            const auto [found, inserted] = by_ends.try_emplace(key, new_id);
            if (inserted) {
                edges.push_back({new_id, {from, to}});
                owners.emplace_back(new_id, a::CanonicalEdge{});
            }
            return found->second;
        };
        for (std::size_t chord = 0; chord < 3; ++chord) {
            for (std::size_t span = 0; span < 4; ++span) {
                const auto index = 4 * chord + span;
                const a::PhysicalFaceId face{100 + index};
                const a::BandId band{50 + span};
                const std::array vertex{point(chord, span), point(chord, span + 1),
                                        point(chord + 1, span + 1), point(chord + 1, span)};
                const std::array outline{edge(vertex[0], vertex[1]), edge(vertex[1], vertex[2]),
                                         edge(vertex[2], vertex[3]), edge(vertex[3], vertex[0])};
                faces.push_back({face, band});
                roles.push_back(
                    {face, chord == 2 ? a::AeroFaceRole::Wake : a::AeroFaceRole::Bound});
                quads.emplace_back(
                    face, vertex, outline, edge(vertex[0], vertex[2]),
                    std::array{a::FacetId{200 + 2 * index}, a::FacetId{201 + 2 * index}});
                if (chord == 2)
                    trailing.push_back({outline[0], {vertex[0], vertex[1]}});
            }
        }
    }

    std::vector<a::PointId> ids;
    std::vector<a::PointBinding> bindings;
    std::vector<double> xyz;
    std::vector<a::EdgeInput> edges;
    std::vector<a::EdgeBindingInput> owners;
    std::vector<a::FaceInput> faces;
    std::vector<a::FaceRoleInput> roles;
    std::vector<a::QuadInput> quads;
    std::vector<a::TrailingEdgeInput> trailing;
};

JsonToken finite_wing_case(const nlohmann::json& wing, double core_radius) {
    auto recipe = root("aerodynamic_steady");
    recipe.field("model").set("geometry", parse(wing.at("case_geometry").dump()));
    recipe.field("model").set("aero_surfaces", parse(wing.at("case_surfaces").dump()));
    recipe.field("flow").set("density", parse("1.2"));
    recipe.field("flow").set(
        "reference",
        parse(
            R"({"density":1.2,"speed":20,"area":8.4,"chord":1.114285714285714,"span":8,"moment_origin":[0.375,0,0]})"));
    recipe.field("flow").field("uniform").set("value", parse(wing.at("flow_m_s").dump()));
    recipe.field("analysis").field("speed").set("value", parse("20"));
    recipe.field("observations").elements[0].set("surfaces", parse(R"(["1","2"])"));
    auto& wake = recipe.field("analysis").field("wake");
    wake.set("core_radius", parse(nlohmann::json(core_radius).dump()));
    wake.set("convection",
             parse(R"({"kind":"prescribed","velocity":{"kind":"constant","value":[20,0,0]}})"));
    return recipe;
}

using Exact = structure::ExactRational;
constexpr std::size_t columns = 8;
using Point = std::array<double, 3>;
using Circulation = std::array<a::ScalarWords, columns>;
using Matrix = std::array<double, columns * columns>;
namespace ls = numerics::linear;

Rat rational_from_view(p::RationalView words) {
    const auto load = [](std::span<const std::uint32_t> limbs) {
        Big result = 0;
        for (const auto word : limbs | std::views::reverse) {
            result <<= 32U;
            result += word;
        }
        return result;
    };
    return Rat{load(words.numerator) * words.sign} / load(words.denominator);
}

Rat as_rational(const Exact& exact, core::ResourcePlan& plan) {
    const auto owner = take(io::detail::copy_rational_words(exact, plan, {}));
    const auto words = owner.view();
    return rational_from_view({words.sign, words.numerator, words.denominator});
}
template <class Field>
Exact exact_quantity(const io::SiCase& si, Field field, core::ResourcePlan& plan) {
    return take(io::materialize_rational(take(si.quantity(field)).recipe, exact_limits(), plan));
}

Exact mix(const Exact& low, const Exact& high, const Exact& alpha, core::ResourcePlan& plan) {
    const auto bound = exact_limits();
    const auto width = take(high.subtract(low, bound, plan));
    return take(low.add(take(width.multiply(alpha, bound, plan)), bound, plan));
}

Point project(const std::array<Exact, 3>& point, core::ResourcePlan& plan) {
    return {as_rational(point[0], plan).convert_to<double>(),
            as_rational(point[1], plan).convert_to<double>(),
            as_rational(point[2], plan).convert_to<double>()};
}

cf::NonpenetrationVector xyz(const Point& point) {
    return {a::SingleWord{point[0]}, a::SingleWord{point[1]}, a::SingleWord{point[2]}};
}

struct CaseLattice {
    std::array<Point, 20> vertices{};
    std::array<Point, columns> collocation{};
    std::array<cf::ExplicitCollocationTangents, columns> tangents{
        cf::ExplicitCollocationTangents{xyz({}), xyz({})},
        cf::ExplicitCollocationTangents{xyz({}), xyz({})},
        cf::ExplicitCollocationTangents{xyz({}), xyz({})},
        cf::ExplicitCollocationTangents{xyz({}), xyz({})},
        cf::ExplicitCollocationTangents{xyz({}), xyz({})},
        cf::ExplicitCollocationTangents{xyz({}), xyz({})},
        cf::ExplicitCollocationTangents{xyz({}), xyz({})},
        cf::ExplicitCollocationTangents{xyz({}), xyz({})}};
};

CaseLattice lattice(const io::SiCase& si, double wake_length, core::ResourcePlan& plan) {
    const auto root_case = take(take(si.source_case()).get().root());
    const auto model = root_case.model.record();
    const structure::PgReferenceLimits bounds{exact_limits(), 1'000'000, 1'000'000, {}};
    const auto make = [&](std::string_view number) {
        return take(Exact::from_decimal(number, exact_limits(), plan));
    };
    const auto center = make("0.5"), quarter = make("0.25"), three_quarters = make("0.75");
    CaseLattice result;
    REQUIRE(model.aero_surfaces.size() == 2);
    REQUIRE(model.geometry.size() == 2);
    for (core::Index side = 0; side < 2; ++side) {
        const auto shape = take(model.aero_surfaces.at(side)).record();
        REQUIRE(std::holds_alternative<io::schema::AeroAttachmentFixed>(shape.attachment.record()));
        REQUIRE(shape.chord.size() == 3);
        REQUIRE(shape.span.size() == 3);
        REQUIRE(shape.sharing.size() == 1);
        REQUIRE(shape.orientation.text() == "positive");
        const auto definition = take(model.geometry.at(side)).record();
        const auto graph = std::get<io::schema::GeometrySurfaceGraph>(definition);
        const auto geometry = shape.geometry.unsigned_value();
        REQUIRE(geometry == graph.id.unsigned_value());
        const auto chart = take(cf::CaseReferenceSurface::prepare(si, geometry, bounds, plan));
        const auto boundary = [&](core::Index id) {
            const auto values = take(graph.boundaries.at(id)).record().vertices;
            REQUIRE(values.size() == 2);
            const auto first = take(values.at(0)), last = take(values.at(1));
            CHECK(as_rational(exact_quantity(si, take(first.at(1)), plan), plan) == 0);
            CHECK(as_rational(exact_quantity(si, take(last.at(1)), plan), plan) == 1);
            return std::array{exact_quantity(si, take(first.at(0)), plan),
                              exact_quantity(si, take(last.at(0)), plan)};
        };
        const auto leading = boundary(static_cast<core::Index>(shape.leading.unsigned_value()));
        const auto trailing = boundary(static_cast<core::Index>(shape.trailing.unsigned_value()));
        const auto sample = [&](const Exact& xi, const Exact& eta) {
            const auto left = mix(leading[0], leading[1], eta, plan);
            const auto right = mix(trailing[0], trailing[1], eta, plan);
            return take(chart.evaluate(mix(left, right, xi, plan), eta, bounds, plan));
        };
        const auto chord = [&](core::Index i) {
            return exact_quantity(si, take(shape.chord.at(i)), plan);
        };
        const auto span = [&](core::Index i) {
            return exact_quantity(si, take(shape.span.at(i)), plan);
        };
        for (core::Index c = 0; c < 3; ++c) {
            const auto first = chord(c);
            const auto xi = c == 2 ? take(first.clone(exact_limits(), plan))
                                   : mix(first, chord(c + 1), quarter, plan);
            for (core::Index s = 0; s < 3; ++s) {
                const auto location = project(sample(xi, span(s)).position, plan);
                const auto index = static_cast<std::size_t>(5 * c + 2 * side + s);
                if (side == 1 && s == 0) {
                    CHECK(location == result.vertices[index]);
                } else {
                    result.vertices[index] = location;
                }
            }
        }
        for (core::Index c = 0; c < 2; ++c) {
            const auto xi = mix(chord(c), chord(c + 1), three_quarters, plan);
            for (core::Index s = 0; s < 2; ++s) {
                const auto point = sample(xi, mix(span(s), span(s + 1), center, plan));
                const auto index = static_cast<std::size_t>(4 * c + 2 * side + s);
                result.collocation[index] = project(point.position, plan);
                result.tangents[index] = {xyz(project(point.chord_tangent, plan)),
                                          xyz(project(point.span_tangent, plan))};
            }
        }
    }
    for (std::size_t i = 0; i < 5; ++i) {
        result.vertices[15 + i] = result.vertices[10 + i];
        result.vertices[15 + i][0] += wake_length;
    }
    return result;
}

ls::Context context(const a::GeometryContext& geometry) {
    const auto trial = geometry.diagnostics.trial;
    if (!trial)
        FAIL("The finite-wing solve requires a trial");
    return {{geometry.fingerprint, geometry.model, geometry.current, core::ChartGeneration{1},
             core::SparsityGeneration{2}, core::NumericalGeneration{3}, core::ScaleGeneration{4},
             geometry.precision, core::Precision::Binary64, core::WakePolicyGeneration{5},
             core::TimeIntervalGeneration{6}, core::BackendGeneration{7}},
            *trial,
            core::ScopeId{8},
            core::ScopeId{9},
            geometry.time};
}

class AlgebraicOnly final : public ls::ConsumerEvaluator<double> {
public:
    [[nodiscard]] numerics::SupportedPrecision supported_precision() const noexcept override {
        return numerics::SupportedPrecision::Binary64AndTwoWord;
    }
    [[nodiscard]] core::Result<core::ErrorCertificate>
    assess(const numerics::PrecisionBuffer<double>&, const numerics::PrecisionBuffer<double>&,
           const ls::Context&, const ls::ConsumerCheck& check,
           core::CertificateContext expected) const override {
        FAIL_CHECK("The explicit original residual follows the dense solve");
        return core::ErrorCertificate::unevaluated(expected, check.budget);
    }
};

a::ScalarWords words(const ls::Solution<double>& solved, std::size_t i) {
    return std::visit(
        [&](const auto& stored) -> a::ScalarWords {
            using T = std::decay_t<decltype(stored)>;
            if constexpr (std::is_same_v<T, numerics::Binary64<double>>)
                return a::SingleWord{stored.values.values()[i]};
            else
                return a::DoubleWord{stored.hi.values()[i], stored.lo.values()[i]};
        },
        solved.values().storage());
}

Rat exact_word(const a::ScalarWords& words) {
    return std::visit(
        [](const auto& value) -> Rat {
            if constexpr (std::is_same_v<std::decay_t<decltype(value)>, a::SingleWord>)
                return rational_word(value.value);
            else
                return rational_word(value.high) + rational_word(value.low);
        },
        words);
}

struct Enclosure {
    Rat center;
    Rat radius;
};

Enclosure unpack(const a::ScaledEnclosure& value) {
    const auto scale = two_power(value.exponent());
    return {exact_word(value.center()) * scale, rational_word(value.radius()) * scale};
}

void check_within(const Enclosure& value, const Rat& reference, const Rat& tolerance) {
    CHECK(abs(value.center - reference) + value.radius <= tolerance);
}

struct SolutionSetup {
    std::array<double, columns> scale{};
    std::array<ls::UnitDimensions, columns> gamma_units{}, speed_units{};
    ls::Context identity;
    ls::Limits solver{8, 1, 0, 5, ls::EstimateAdmission{1U << 20U}};
    ls::Accuracy accuracy{ls::Norm::Infinity,     1e-25, 0., 1e-23, core::ObservationId{20},
                          core::ObservationId{21}};

    explicit SolutionSetup(const a::GeometryContext& geometry) : identity(context(geometry)) {
        scale.fill(20.);
        gamma_units.fill({2, 0, -1, 0, 0, 0, 0});
        speed_units.fill({1, 0, -1, 0, 0, 0, 0});
    }
    [[nodiscard]] ls::Scales scales() const {
        return {scale, scale, gamma_units, speed_units, identity};
    }
};

Circulation solve(const Matrix& values, const std::array<double, columns>& rhs,
                  const SolutionSetup& setup, core::ResourcePlan& plan) {
    constexpr auto operation = core::Operation::SolveLinear;
    const auto matrix =
        take(numerics::DenseMatrix<double>::copy({8, 8, 8}, values, plan, operation));
    auto source = take(numerics::DenseMatrix<double>::copy({8, 1, 8}, rhs, plan, operation));
    const auto right = take(numerics::PrecisionBuffer<double>::binary64(
        std::move(source), core::NumericalGeneration{30}));
    const ls::DenseOperator op{matrix, setup.identity};
    const auto scales = setup.scales();
    auto prepared = take(ls::DensePrepared::prepare(op, scales, setup.solver, plan));
    const AlgebraicOnly algebraic;
    const ls::Precision<double> precision{
        ls::RetainedPrecision{core::NumericalGeneration{31}}, {}, algebraic};
    const auto solved =
        take(ls::solve_linear(op, ls::Rhs{right, setup.identity, core::RhsGeneration{40}}, scales,
                              precision, setup.accuracy, setup.solver, prepared, plan));
    Circulation result{a::SingleWord{0}, a::SingleWord{0}, a::SingleWord{0}, a::SingleWord{0},
                       a::SingleWord{0}, a::SingleWord{0}, a::SingleWord{0}, a::SingleWord{0}};
    for (std::size_t i = 0; i < columns; ++i)
        result[i] = words(solved, i);
    return result;
}

constexpr std::size_t band_count = 20;
using EnclosedVector = std::array<Enclosure, 3>;
using Bands = std::array<EnclosedVector, band_count>;

Enclosure add(const Enclosure& a, const Enclosure& b) {
    return {a.center + b.center, a.radius + b.radius};
}

Enclosure subtract(const Enclosure& a, const Enclosure& b) {
    return {a.center - b.center, a.radius + b.radius};
}

Enclosure multiply(const Enclosure& a, const Enclosure& b) {
    return {a.center * b.center,
            abs(a.center) * b.radius + abs(b.center) * a.radius + a.radius * b.radius};
}

Enclosure scale(const Enclosure& value, const Rat& factor) {
    return {value.center * factor, value.radius * abs(factor)};
}

EnclosedVector unpack_vector(const cf::DiscreteLoadVector& value) {
    return {unpack(value[0]), unpack(value[1]), unpack(value[2])};
}

EnclosedVector arm_moment(EnclosedVector point, const Point& origin, const EnclosedVector& force) {
    for (std::size_t i = 0; i < point.size(); ++i)
        point[i].center -= rational_word(origin[i]);
    EnclosedVector result;
    for (std::size_t i = 0; i < result.size(); ++i) {
        const auto j = (i + 1) % result.size(), k = (i + 2) % result.size();
        result[i] = subtract(multiply(point[j], force[k]), multiply(point[k], force[j]));
    }
    return result;
}

using RationalPoint = std::array<Rat, 3>;

struct WrenchBands {
    Bands force, moment;
};

Rat band_border(std::size_t i) {
    return Rat{-4} + Rat{2} * i / 5;
}

RationalPoint exact_position(const Point& words) {
    return {rational_word(words[0]), rational_word(words[1]), rational_word(words[2])};
}

void accumulate(WrenchBands& bands, std::size_t band, const EnclosedVector& force,
                const EnclosedVector& moment, const Rat& weight) {
    for (std::size_t i = 0; i < force.size(); ++i) {
        bands.force[band][i] = add(bands.force[band][i], scale(force[i], weight));
        bands.moment[band][i] = add(bands.moment[band][i], scale(moment[i], weight));
    }
}

void integrate_edge(WrenchBands& bands, const RationalPoint& first, const RationalPoint& last,
                    const EnclosedVector& point, const EnclosedVector& force, const Point& origin) {
    REQUIRE(first[1] >= band_border(0));
    REQUIRE(last[1] >= band_border(0));
    REQUIRE(first[1] <= band_border(band_count));
    REQUIRE(last[1] <= band_border(band_count));
    if (first[1] == last[1]) {
        // The flow-aligned bound edges keep the point-force observation.
        REQUIRE(point[1].radius == 0);
        REQUIRE(point[1].center == first[1]);
        std::size_t band = 0;
        while (band + 1 < band_count && point[1].center >= band_border(band + 1))
            ++band;
        accumulate(bands, band, force, arm_moment(point, origin, force), Rat{1});
        return;
    }
    RationalPoint delta{};
    auto displaced = point;
    for (std::size_t i = 0; i < delta.size(); ++i) {
        delta[i] = last[i] - first[i];
        displaced[i].center -= (first[i] + last[i]) / 2;
    }
    const auto couple = arm_moment(displaced, Point{}, force);
    Rat total_weight{};
    for (std::size_t band = 0; band < band_count; ++band) {
        const Rat at_low = (band_border(band) - first[1]) / delta[1];
        const Rat at_high = (band_border(band + 1) - first[1]) / delta[1];
        const Rat low = std::max(Rat{0}, std::min(at_low, at_high));
        const Rat high = std::min(Rat{1}, std::max(at_low, at_high));
        if (low >= high)
            continue;
        const Rat middle = (low + high) / 2;
        EnclosedVector location;
        for (std::size_t i = 0; i < location.size(); ++i)
            location[i] = {first[i] + middle * delta[i], Rat{0}};
        auto moment = arm_moment(location, origin, force);
        for (std::size_t i = 0; i < moment.size(); ++i)
            moment[i] = add(moment[i], couple[i]);
        const Rat weight = high - low;
        total_weight += weight;
        accumulate(bands, band, force, moment, weight);
    }
    CHECK(total_weight == 1);
}

[[maybe_unused]] void check_line_invariants() {
    const RationalPoint first{0, -2, 0}, last{4, 2, 0};
    const EnclosedVector load_point{{{1, 0}, {-1, 0}, {0, 0}}};
    const EnclosedVector force{{{0, 0}, {0, 0}, {8, 0}}};
    WrenchBands forward{}, backward{}, shifted{};
    integrate_edge(forward, first, last, load_point, force, Point{});
    integrate_edge(backward, last, first, load_point, force, Point{});
    const Point shift{1., -.5, .25};
    integrate_edge(shifted, first, last, load_point, force, shift);
    EnclosedVector total_force{}, total_moment{};
    for (std::size_t band = 0; band < band_count; ++band) {
        CAPTURE(band);
        const EnclosedVector offset{{{1, 0}, {-Rat{1} / 2, 0}, {Rat{1} / 4, 0}}};
        const auto transport = arm_moment(offset, Point{}, forward.force[band]);
        for (std::size_t i = 0; i < force.size(); ++i) {
            check_within(subtract(forward.force[band][i], backward.force[band][i]), 0, 0);
            check_within(subtract(forward.moment[band][i], backward.moment[band][i]), 0, 0);
            check_within(
                subtract(shifted.moment[band][i], subtract(forward.moment[band][i], transport[i])),
                0, 0);
            total_force[i] = add(total_force[i], forward.force[band][i]);
            total_moment[i] = add(total_moment[i], forward.moment[band][i]);
        }
    }
    const auto point_moment = arm_moment(load_point, Point{}, force);
    for (std::size_t i = 0; i < force.size(); ++i) {
        check_within(subtract(total_force[i], force[i]), 0, 0);
        check_within(subtract(total_moment[i], point_moment[i]), 0, 0);
    }
    // The compensating couple shifts the moment in a strip without moving its force.
    check_within(forward.force[7][2], Rat{4} / 5, 0);
    check_within(forward.moment[7][0], -Rat{8} / 5, 0);
    check_within(forward.moment[7][1], 0, 0);
}

RationalPoint endpoint(std::span<const Point> vertices, a::PointId id) {
    REQUIRE(id.value() > 0);
    REQUIRE(id.value() <= vertices.size());
    return exact_position(vertices[static_cast<std::size_t>(id.value() - 1)]);
}

std::array<Enclosure, band_count>
check_span_lines(const cf::DiscreteCirculationLoads& loads, std::span<const Point> vertices,
                 const nlohmann::json& expected, const Point& origin, const Rat& cos_alpha,
                 const Rat& sin_alpha, const Rat& force_tolerance, const Rat& moment_tolerance) {
    std::array<bool, 18> found{};
    const auto& edges = expected.at("edges");
    REQUIRE(loads.edge_count() == found.size());
    REQUIRE(edges.size() == found.size());
    WrenchBands observed{};
    for (std::size_t edge = 0; edge < loads.edge_count(); ++edge) {
        const auto load = take(loads.edge(edge));
        CAPTURE(load.edge.value());
        auto pair = std::array{load.endpoints[0].value(), load.endpoints[1].value()};
        std::ranges::sort(pair);
        const auto reference = std::ranges::find_if(edges, [&](const auto& item) {
            return item.at("endpoint_ids").template get<decltype(pair)>() == pair;
        });
        if (reference == edges.end())
            FAIL("Every original load edge must have one independent reference");
        const auto id = static_cast<std::size_t>(reference - edges.begin());
        REQUIRE_FALSE(found.at(id));
        found.at(id) = true;
        const auto point = unpack_vector(load.position);
        const auto force = unpack_vector(load.force);
        const auto first = endpoint(vertices, load.endpoints[0]);
        const auto last = endpoint(vertices, load.endpoints[1]);
        CHECK(load.excluded_bound_occurrences > 0);
        for (std::size_t i = 0; i < force.size(); ++i) {
            const auto oracle_point = rational_word(reference->at("point_m").at(i).get<double>());
            const auto oracle_force = rational_word(reference->at("force_N").at(i).get<double>());
            check_within(point[i], oracle_point, rational_word(2e-14));
            check_within(point[i], (first[i] + last[i]) / 2, rational_word(2e-14));
            check_within(force[i], oracle_force, force_tolerance);
        }
        integrate_edge(observed, first, last, point, force, origin);
    }
    std::array<Enclosure, band_count> lift{};
    EnclosedVector sum_force{}, sum_moment{};
    for (std::size_t band = 0; band < band_count; ++band) {
        CAPTURE(band);
        lift[band] = add(scale(observed.force[band][0], -sin_alpha),
                         scale(observed.force[band][2], cos_alpha));
        check_within(lift[band], rational_word(expected.at("span_lift_N").at(band).get<double>()),
                     force_tolerance);
        for (std::size_t i = 0; i < sum_force.size(); ++i) {
            CAPTURE(i);
            const auto f = rational_word(expected.at("span_forces_N").at(band).at(i).get<double>());
            const auto m =
                rational_word(expected.at("span_moments_Nm").at(band).at(i).get<double>());
            check_within(observed.force[band][i], f, force_tolerance);
            check_within(observed.moment[band][i], m, moment_tolerance);
            sum_force[i] = add(sum_force[i], observed.force[band][i]);
            sum_moment[i] = add(sum_moment[i], observed.moment[band][i]);
        }
    }
    const auto total_force = unpack_vector(loads.total_force());
    const auto total_moment = unpack_vector(loads.total_moment());
    for (std::size_t i = 0; i < sum_force.size(); ++i) {
        check_within(subtract(sum_force[i], total_force[i]), 0, rational_word(1e-10));
        check_within(subtract(sum_moment[i], total_moment[i]), 0, rational_word(1e-10));
    }
    return lift;
}

Rat trefftz_midpoint(const cf::DiscreteCirculation& solved, const FiniteWingMesh& mesh,
                     const a::GlobalContactState& contact, std::span<const Point> vertices,
                     const Rat& density) {
    constexpr std::size_t strips = 4;
    std::array<RationalPoint, strips + 1> edge{};
    std::array<Rat, strips> trace{};
    for (std::size_t i = 0; i < strips; ++i) {
        const auto wake =
            take_global(contact.uvlm_quad(mesh.quads[columns + i].facets[0], 1000000));
        const auto& owner = wake.attachment;
        if (!owner)
            FAIL("The Trefftz trace must come from a qualified trailing owner");
        const auto first = endpoint(vertices, wake.corners[0]);
        if (i != 0)
            CHECK(edge[i] == first);
        edge[i] = first;
        edge[i + 1] = endpoint(vertices, wake.corners[1]);
        const auto row = columns - strips + i;
        REQUIRE(take(solved.row(row)).face == owner->bound_face);
        trace[i] = exact_word(words(solved.solution(), row)) * owner->coefficient;
    }
    std::array<Rat, strips + 1> jumps{};
    for (std::size_t i = 0; i < strips; ++i) {
        jumps[i] -= trace[i];
        jumps[i + 1] += trace[i];
    }
    Rat sum{};
    for (std::size_t i = 0; i < strips; ++i) {
        const Rat dy = edge[i + 1][1] - edge[i][1];
        const Rat dz = edge[i + 1][2] - edge[i][2];
        REQUIRE((dy != 0 || dz != 0));
        const Rat y = edge[i][1] + dy / 2;
        const Rat z = edge[i][2] + dz / 2;
        for (std::size_t j = 0; j < edge.size(); ++j) {
            const Rat ry = y - edge[j][1], rz = z - edge[j][2];
            const Rat r2 = ry * ry + rz * rz;
            REQUIRE(r2 > 0);
            sum += trace[i] * jumps[j] * (ry * dy + rz * dz) / r2;
        }
    }
    return -density * sum / (4 * rational_word(std::numbers::pi));
}

struct SpanObservation {
    Enclosure drag;
    Enclosure lift;
    Enclosure pitching_moment;
    std::array<Enclosure, band_count> band_lift;
};

SpanObservation observe(const cf::DiscreteCirculationLoads& loads, const nlohmann::json& expected,
                        const Point& origin, const Rat& cos_alpha, const Rat& sin_alpha,
                        const Rat& force_tolerance, const Rat& moment_tolerance) {
    const auto lift = [&](const EnclosedVector& value) {
        return add(scale(value[0], -sin_alpha), scale(value[2], cos_alpha));
    };
    const auto drag = [&](const EnclosedVector& value) {
        return add(scale(value[0], cos_alpha), scale(value[2], sin_alpha));
    };
    REQUIRE(loads.edge_count() == 18);
    Bands force_bands{}, moment_bands{};
    for (std::size_t edge = 0; edge < loads.edge_count(); ++edge) {
        const auto value = take(loads.edge(edge));
        const auto point = unpack_vector(value.position);
        const auto force = unpack_vector(value.force);
        const auto moment = arm_moment(point, origin, force);
        const auto& y = point[1];
        // These frozen midpoints have exact dyadic y words. Never classify a
        // range crossing a strip boundary by rounding its centre.
        REQUIRE(y.radius == 0);
        REQUIRE(y.center >= -4);
        REQUIRE(y.center <= 4);
        std::size_t band = 0;
        while (band + 1 < band_count && y.center >= band_border(band + 1))
            ++band;
        for (std::size_t k = 0; k < force.size(); ++k) {
            force_bands[band][k] = add(force_bands[band][k], force[k]);
            moment_bands[band][k] = add(moment_bands[band][k], moment[k]);
        }
    }
    const auto force = unpack_vector(loads.total_force());
    const auto moment = unpack_vector(loads.total_moment());
    EnclosedVector force_sum{}, moment_sum{};
    SpanObservation result{drag(force), lift(force), moment[1], {}};
    check_within(result.drag, rational_word(expected.at("drag_lift_N").at(0).get<double>()),
                 force_tolerance);
    check_within(result.lift, rational_word(expected.at("drag_lift_N").at(1).get<double>()),
                 force_tolerance);
    for (std::size_t band = 0; band < band_count; ++band) {
        CAPTURE(band);
        result.band_lift[band] = lift(force_bands[band]);
        check_within(result.band_lift[band],
                     rational_word(expected.at("band_lift_N").at(band).get<double>()),
                     force_tolerance);
        for (std::size_t k = 0; k < force.size(); ++k) {
            CAPTURE(k);
            force_sum[k] = add(force_sum[k], force_bands[band][k]);
            moment_sum[k] = add(moment_sum[k], moment_bands[band][k]);
            check_within(force_bands[band][k],
                         rational_word(expected.at("band_forces_N").at(band).at(k).get<double>()),
                         force_tolerance);
            check_within(moment_bands[band][k],
                         rational_word(expected.at("band_moments_Nm").at(band).at(k).get<double>()),
                         moment_tolerance);
        }
    }
    const auto conservation = rational_word(1e-10);
    for (std::size_t axis = 0; axis < force.size(); ++axis) {
        CAPTURE(axis);
        check_within(force[axis], rational_word(expected.at("force_N").at(axis).get<double>()),
                     force_tolerance);
        check_within(moment[axis], rational_word(expected.at("moment_Nm").at(axis).get<double>()),
                     moment_tolerance);
        check_within(subtract(force_sum[axis], force[axis]), Rat{0}, conservation);
        check_within(subtract(moment_sum[axis], moment[axis]), Rat{0}, conservation);
    }
    return result;
}

void compare_ideal(const SpanObservation& actual, const nlohmann::json& finite,
                   const nlohmann::json& ideal, const Rat& force_tolerance,
                   const Rat& moment_tolerance) {
    const auto compare = [&](const Enclosure& computed, double finite_value, double ideal_value,
                             const Rat& tolerance) {
        const Rat exact_ideal = rational_word(ideal_value);
        const Rat estimate = rational_word(finite_value) - exact_ideal;
        check_within(computed, exact_ideal, abs(estimate) + tolerance);
        if (abs(estimate) > tolerance) {
            const auto difference =
                scale(subtract(computed, {exact_ideal, Rat{0}}), estimate > 0 ? Rat{1} : Rat{-1});
            CHECK(difference.center - difference.radius > 0);
        }
    };
    compare(actual.drag, finite.at("drag_lift_N").at(0).get<double>(),
            ideal.at("drag_lift_N").at(0).get<double>(), force_tolerance);
    compare(actual.lift, finite.at("drag_lift_N").at(1).get<double>(),
            ideal.at("drag_lift_N").at(1).get<double>(), force_tolerance);
    compare(actual.pitching_moment, finite.at("moment_Nm").at(1).get<double>(),
            ideal.at("moment_Nm").at(1).get<double>(), moment_tolerance);
    for (std::size_t band = 0; band < band_count; ++band) {
        CAPTURE(band);
        compare(actual.band_lift[band], finite.at("band_lift_N").at(band).get<double>(),
                ideal.at("band_lift_N").at(band).get<double>(), force_tolerance);
    }
}

struct WingCondition {
    double wake_length{};
    double core_radius{};
    bool full_polar{};
    const nlohmann::json* ideal_angle{};
    const nlohmann::json* span_reference{};
};

void run_wing(const nlohmann::json& wing, const nlohmann::json& expected,
              const nlohmann::json& reference_metadata, const nlohmann::json& polar,
              const WingCondition& condition) {
    auto recipe = finite_wing_case(wing, condition.core_radius);
    const auto si = materialize(recipe, io::SourceId{532});
    core::ResourcePlan memory{1024ULL * 1024 * 1024};
    const auto case_grid = lattice(si, condition.wake_length, memory);
    for (std::size_t i = 0; i < case_grid.vertices.size(); ++i) {
        const bool translated = !condition.full_polar && i >= 15;
        auto reference = wing.at("vertices_m").at(translated ? i - 5 : i).get<Point>();
        if (translated)
            reference[0] += condition.wake_length;
        for (std::size_t axis = 0; axis < 3; ++axis)
            CHECK(std::abs(case_grid.vertices[i][axis] - reference[axis]) < 2e-14);
    }
    for (std::size_t i = 0; i < case_grid.collocation.size(); ++i)
        for (std::size_t axis = 0; axis < 3; ++axis)
            CHECK(std::abs(case_grid.collocation[i][axis] -
                           wing.at("collocations_m").at(i).at(axis).get<double>()) < 2e-14);
    const FiniteWingMesh mesh{case_grid.vertices};
    const std::array bands{a::BandId{50}, a::BandId{51}, a::BandId{52}, a::BandId{53}};
    auto geometry_caps = test_geometry::limits();
    geometry_caps.max_steps = 1000000;
    auto geometry_work = test_geometry::workspace(memory);
    const Number instant{5, 4};
    const auto typed = take(take(si.source_case()).get().root());
    const auto analysis = std::get<io::schema::AnalysisStatic>(typed.analysis.record());
    const auto wake_input = std::get<io::schema::WakeActive>(analysis.wake.record());
    const auto radius = as_rational(exact_quantity(si, wake_input.core_radius, memory), memory)
                            .convert_to<double>();
    CHECK(radius == condition.core_radius);
    const a::SinglePositions positions{mesh.xyz};
    const a::GeometryInput input_geometry{mesh.ids,
                                          mesh.bindings,
                                          positions,
                                          positions,
                                          mesh.edges,
                                          mesh.faces,
                                          bands,
                                          {},
                                          mesh.quads,
                                          a::SingleWord{radius},
                                          test_geometry::context(),
                                          instant.view()};
    const auto geometry = test_geometry::take(
        a::prepare_owned_geometry(input_geometry, geometry_caps, *geometry_work, memory));
    const a::GlobalLimits contact_caps{std::numeric_limits<std::uint64_t>::max() / 4,
                                       ls::EstimateAdmission{65536}};
    auto global = take(a::GlobalWorkspace::create(memory, contact_caps.metadata));
    const auto topology = take_global(a::build_global_topology(
        geometry, {mesh.roles, mesh.owners, mesh.trailing}, contact_caps, *global, memory));
    const auto contact =
        take_global(a::check_global_geometry(topology, geometry, contact_caps, *global, memory));
    REQUIRE(topology.trailing_edges().size() == 4);

    std::vector<a::ResolvedQuery> queries;
    queries.reserve(columns);
    for (std::size_t i = 0; i < columns; ++i)
        queries.push_back(test_geometry::take(a::resolve_query(
            geometry, {a::QueryId{i + 1}, a::FreePoint{a::FreeSingle{case_grid.collocation[i]}}},
            geometry.context(), geometry_caps, *geometry_work, memory)));

    const cf::DiscreteNonpenetrationLimits caps{
        {core::Precision::TwoWord, 1e-12},
        {core::Precision::TwoWord, 1e-12, 1e-12},
        {geometry_caps, {ls::EstimateAdmission{65536}, 1000000000}, UINT64_C(32000000000)},
        1000000,
        UINT64_C(1000000000000)};
    a::UniformFlowWorkspace flow_work;
    const auto flow = take(cf::prepare_case_uniform_flow(si, stamp, {}, flow_work, limits, memory));
    const auto projected = [&](const cf::CaseUniformFlowReady& field) {
        auto sampled = take(cf::evaluate_case_uniform_flow(field, field, {stamp, instant.view()},
                                                           flow_work, limits.flow, memory));
        const Number one{1}, absolute{1, 1000000000}, zero{0};
        const cf::CaseFlowComponentBudget component{one.view(), absolute.view(), zero.view()};
        const std::array bounds{component, component, component};
        return take(cf::project_case_uniform_flow(
            std::move(sampled), {core::Precision::TwoWord, bounds, bounds, {4096, budget}}));
    };

    const auto parameters = [&](const Circulation& gamma) {
        std::vector<cf::DiscreteCollocation> rows;
        rows.reserve(columns);
        for (std::size_t i = 0; i < columns; ++i)
            rows.emplace_back(mesh.quads[i].facets[0], queries[i], case_grid.tangents[i], xyz({}),
                              gamma[i], a::SingleWord{0});
        std::vector<cf::FrozenWakeTrace> traces;
        std::vector<cf::DiscreteWakePanel> wake;
        for (std::size_t i = 0; i < 4; ++i) {
            const a::P1StrengthId history{700 + i};
            traces.push_back({history, gamma[4 + i]});
            wake.emplace_back(mesh.quads[8 + i].facets[0], cf::AttachedWakeFront{}, history);
        }
        return std::tuple{std::move(rows), std::move(traces), std::move(wake)};
    };
    constexpr a::P1ActionStamp trial{815, 6};
    a::SharedP1Workspace shared;
    auto local = take(a::LocalPacketWorkspace::create(memory));
    auto kernel = take(a::P1Workspace::create(memory));
    const auto evaluate = [&](const Circulation& values) {
        const auto [rows, traces, wake] = parameters(values);
        const cf::DiscreteNonpenetrationInput input{contact,
                                                    rows,
                                                    wake,
                                                    {core::AcceptedStateId{2}, traces},
                                                    trial,
                                                    stamp,
                                                    p::RequestedSide::Automatic};
        return take(cf::evaluate_discrete_nonpenetration(input, projected(flow), caps, shared,
                                                         *local, *kernel, memory));
    };
    Circulation basis{a::SingleWord{0}, a::SingleWord{0}, a::SingleWord{0}, a::SingleWord{0},
                      a::SingleWord{0}, a::SingleWord{0}, a::SingleWord{0}, a::SingleWord{0}};
    const auto zero = evaluate(basis);
    std::array<Enclosure, columns> offsets;
    std::array<double, columns> rhs{};
    const auto tolerance = rational_word(2e-10);
    for (std::size_t i = 0; i < columns; ++i) {
        const auto row = take(zero.row(i));
        CHECK(row.face == mesh.faces[i].id);
        offsets[i] = unpack(row.residual);
        rhs[i] = (-offsets[i].center).convert_to<double>();
        check_within({-offsets[i].center, offsets[i].radius},
                     rational_word(expected.at("rhs_m_s").at(i).get<double>()), tolerance);
    }
    Matrix matrix{};
    for (std::size_t j = 0; j < columns; ++j) {
        CAPTURE(j);
        basis[j] = a::SingleWord{1};
        const auto sampled = evaluate(basis);
        basis[j] = a::SingleWord{0};
        for (std::size_t i = 0; i < columns; ++i) {
            CAPTURE(i);
            const auto coefficient = unpack(take(sampled.row(i)).residual);
            const Enclosure assembled{coefficient.center - offsets[i].center,
                                      coefficient.radius + offsets[i].radius};
            matrix[j * columns + i] = assembled.center.convert_to<double>();
            check_within(
                assembled,
                rational_word(expected.at("steady_columns_per_m").at(i).at(j).get<double>()),
                tolerance);
        }
    }
    const SolutionSetup setup{geometry.context()};
    const auto base_gamma = solve(matrix, rhs, setup, memory);
    const auto residual = evaluate(base_gamma);
    for (std::size_t i = 0; i < columns; ++i) {
        CAPTURE(i);
        check_within(unpack(take(residual.row(i)).residual), Rat{0}, tolerance);
        CHECK(abs(exact_word(base_gamma[i]) -
                  rational_word(expected.at("circulation_m2_s").at(i).get<double>())) <= tolerance);
    }
    // Release the same steady circulation into the fixed downstream trace, then
    // verify the loads path's own circulation solve and native original residual.
    std::vector<cf::BoundCirculationRate> rates;
    std::vector<a::EdgeId> physical_edges;
    for (std::size_t i = 0; i < columns; ++i) {
        rates.emplace_back(mesh.faces[i].id, a::SingleWord{0});
        const auto& edges = mesh.quads[i].edges;
        physical_edges.insert(physical_edges.end(), edges.begin(), edges.end());
    }
    std::ranges::sort(physical_edges);
    const auto duplicate = std::ranges::unique(physical_edges);
    physical_edges.erase(duplicate.begin(), duplicate.end());
    std::vector<cf::EdgeMidpointVelocity> midpoint;
    for (const auto edge : physical_edges) {
        const auto is_trailing = std::ranges::any_of(
            mesh.trailing, [edge](const auto& item) { return item.edge == edge; });
        if (!is_trailing)
            midpoint.emplace_back(edge, xyz({}));
    }
    const auto fluid = std::get<io::schema::FlowKnown>(typed.flow.record());
    const auto reference = fluid.reference.record();
    const auto read_value = [&](auto field) {
        return as_rational(exact_quantity(si, field, memory), memory).template convert_to<double>();
    };
    Point origin{};
    for (std::size_t k = 0; k < origin.size(); ++k)
        origin[k] = read_value(take(reference.moment_origin.at(static_cast<core::Index>(k))));
    const double density = read_value(fluid.density);
    CHECK(origin == reference_metadata.at("moment_origin_m").get<Point>());
    CHECK(density == reference_metadata.at("density_kg_m3").get<double>());
    const cf::DiscreteCirculationRequest request{setup.scales(),
                                                 setup.solver,
                                                 setup.accuracy,
                                                 core::RhsGeneration{3},
                                                 core::NumericalGeneration{4},
                                                 1e-13,
                                                 1e-13,
                                                 1e-11,
                                                 0.};
    const cf::DiscreteVectorLoadInput load_input{trial, request.retained_generation, midpoint,
                                                 rates, a::SingleWord{density},      xyz(origin)};
    const auto& angles = polar.at("angles");
    REQUIRE(angles.size() == (condition.full_polar ? 5U : 1U));
    std::array<SpanObservation, 5> observed{};
    std::array<std::array<Enclosure, band_count>, 5> line_lift{};
    std::array<Rat, 5> far_drag{}, sine{}, cosine{};
    const auto velocity_scale = rational_word(read_value(reference.speed));
    const auto dynamic_area = rational_word(density) * velocity_scale * velocity_scale *
                              rational_word(read_value(reference.area)) / 2;
    const auto chord = rational_word(read_value(reference.chord));
    const auto force_tolerance = rational_word(2e-11) * dynamic_area;
    const auto moment_tolerance = force_tolerance * chord;
    const auto base_projection = projected(flow);
    const auto base_z = rational_from_view(take(base_projection.source().jet.velocity(2)));
    REQUIRE(base_z != 0);

    for (std::size_t i = 0; i < angles.size(); ++i) {
        const auto& sample = angles.at(i);
        const auto alpha = sample.at("alpha_degrees").get<int>();
        CAPTURE(alpha);
        REQUIRE(alpha == (condition.full_polar ? static_cast<int>(i) - 2 : 1));
        auto angled = recipe;
        angled.field("flow").field("uniform").set("value", parse(sample.at("flow_m_s").dump()));
        const auto angled_si = materialize(angled, io::SourceId{532});
        const auto angled_flow =
            take(cf::prepare_case_uniform_flow(angled_si, stamp, {}, flow_work, limits, memory));
        auto ambient = projected(angled_flow);
        const auto ux = rational_from_view(take(ambient.source().jet.velocity(0)));
        const auto uz = rational_from_view(take(ambient.source().jet.velocity(2)));

        // The frozen sections have a chord tangent parallel to +x, hence their
        // normals all have an exact zero x-component. Only U_z changes the RHS.
        std::array<double, columns> angled_rhs{};
        for (std::size_t row = 0; row < columns; ++row) {
            const Rat value = -offsets[row].center * uz / base_z;
            angled_rhs[row] = value.convert_to<double>();
        }
        const auto gamma = solve(matrix, angled_rhs, setup, memory);
        const auto [rows, traces, wake] = parameters(gamma);
        const cf::DiscreteNonpenetrationInput input{contact,
                                                    rows,
                                                    wake,
                                                    {core::AcceptedStateId{2}, traces},
                                                    trial,
                                                    stamp,
                                                    p::RequestedSide::Automatic};
        const auto loads = take(cf::solve_discrete_circulation_loads(
            input, load_input, std::move(ambient), caps, request, *geometry_work, shared, *local,
            *kernel, memory));
        for (std::size_t row = 0; row < columns; ++row) {
            CAPTURE(row);
            CHECK(abs(exact_word(gamma[row]) -
                      rational_word(sample.at("circulation_m2_s").at(row).get<double>())) <=
                  tolerance);
            CHECK(abs(exact_word(words(loads.primal().solution(), row)) - exact_word(gamma[row])) <=
                  tolerance);
            check_within(unpack(take(loads.primal().row(row)).original_residual), Rat{0},
                         tolerance);
        }
        observed[i] = observe(loads, sample, origin, ux / velocity_scale, uz / velocity_scale,
                              force_tolerance, moment_tolerance);
        if (condition.span_reference != nullptr) {
            sine[i] = uz / velocity_scale;
            cosine[i] = ux / velocity_scale;
            const auto& line = condition.span_reference->at("angles").at(i);
            REQUIRE(line.at("alpha_degrees").get<int>() == alpha);
            line_lift[i] =
                check_span_lines(loads, case_grid.vertices, line, origin, ux / velocity_scale,
                                 uz / velocity_scale, force_tolerance, moment_tolerance);
            far_drag[i] = trefftz_midpoint(loads.primal(), mesh, contact, case_grid.vertices,
                                           rational_word(density));
            check_within({far_drag[i], Rat{0}},
                         rational_word(line.at("trefftz_drag_N").get<double>()), force_tolerance);
        }
        if (condition.ideal_angle != nullptr)
            compare_ideal(observed[i], sample, *condition.ideal_angle, force_tolerance,
                          moment_tolerance);
        if (alpha == 0) {
            const auto force = unpack_vector(loads.total_force());
            const auto moment = unpack_vector(loads.total_moment());
            for (std::size_t k = 0; k < force.size(); ++k) {
                check_within(force[k], Rat{0}, rational_word(1e-9) * dynamic_area);
                check_within(moment[k], Rat{0}, rational_word(1e-9) * dynamic_area * chord);
            }
        }
    }

    if (!condition.full_polar)
        return;
    for (std::size_t delta = 1; delta <= 2; ++delta) {
        CAPTURE(delta);
        const auto& sample = polar.at("slopes").at(delta - 1);
        REQUIRE(sample.at("delta_degrees").get<std::size_t>() == delta);
        const auto& minus = observed[2 - delta];
        const auto& plus = observed[2 + delta];
        const Rat denominator = 2 * Rat{delta} * rational_word(std::numbers::pi) / 180;
        const auto coefficient_scale = Rat{1} / (dynamic_area * denominator);
        const auto slope_tolerance = rational_word(4e-11) / denominator;
        const auto slope = [&](const Enclosure& neg, const Enclosure& pos) {
            return scale(subtract(pos, neg), coefficient_scale);
        };
        const auto lift_slope = slope(minus.lift, plus.lift);
        const auto moment_slope =
            scale(slope(minus.pitching_moment, plus.pitching_moment), Rat{1} / chord);
        const auto lift_reference = rational_word(sample.at("cl_per_rad").get<double>());
        const auto moment_reference = rational_word(sample.at("cmy_per_rad").get<double>());
        check_within(lift_slope, lift_reference, slope_tolerance);
        check_within(moment_slope, moment_reference, slope_tolerance);
        check_within(lift_slope, lift_reference, abs(lift_reference) / 100);
        const Rat moment_scale = std::max(Rat{abs(moment_reference)}, Rat{1} / 10);
        check_within(moment_slope, moment_reference, moment_scale / 100);
        Rat norm_squared{};
        Rat error_squared{};
        const Rat inverse_width = Rat{5} / 2;
        for (std::size_t band = 0; band < band_count; ++band) {
            CAPTURE(band);
            const auto band_slope = slope(minus.band_lift[band], plus.band_lift[band]);
            const auto band_reference =
                rational_word(sample.at("band_cl_per_rad").at(band).get<double>());
            const Rat error = abs(band_slope.center - band_reference) + band_slope.radius;
            error_squared += error * error * inverse_width;
            norm_squared += band_reference * band_reference * inverse_width;
            check_within(band_slope, band_reference, slope_tolerance);
        }
        REQUIRE(norm_squared > 0);
        CHECK(error_squared <= norm_squared / 2500);
        if (condition.span_reference != nullptr) {
            const auto& line = condition.span_reference->at("slopes").at(delta - 1);
            REQUIRE(line.at("delta_degrees").get<std::size_t>() == delta);
            Rat reference_norm{}, error_norm{};
            for (std::size_t band = 0; band < band_count; ++band) {
                CAPTURE(band);
                const auto measured = slope(line_lift[2 - delta][band], line_lift[2 + delta][band]);
                const auto band_reference =
                    rational_word(line.at("span_lift_per_rad").at(band).get<double>());
                const Rat error = abs(measured.center - band_reference) + measured.radius;
                error_norm += inverse_width * error * error;
                reference_norm += inverse_width * band_reference * band_reference;
                check_within(measured, band_reference, slope_tolerance);
            }
            REQUIRE(reference_norm > 0);
            CHECK(error_norm <= reference_norm / 2500);
            const Rat radians = denominator / 2;
            const Rat drag_scale = Rat{1} / (2 * radians * radians);
            const Enclosure near =
                scale(subtract(add(minus.drag, plus.drag), scale(observed[2].drag, 2)), drag_scale);
            const Rat far =
                (far_drag[2 - delta] + far_drag[2 + delta] - 2 * far_drag[2]) * drag_scale;
            const Rat drag_tolerance = 4 * force_tolerance * drag_scale;
            check_within(near, rational_word(line.at("near_drag_N_per_rad2").get<double>()),
                         drag_tolerance);
            check_within({far, Rat{0}},
                         rational_word(line.at("trefftz_drag_N_per_rad2").get<double>()),
                         drag_tolerance);
            // Fixed geometry gives Gamma(alpha)=sin(alpha)*Gamma'(0). The
            // ambient term is orthogonal to the wind-axis direction, and
            // the opposite samples cancel the induced z-force projection.
            REQUIRE(sine[2 - delta] == -sine[2 + delta]);
            REQUIRE(cosine[2 - delta] == cosine[2 + delta]);
            const Rat symmetric_sine = 2 * sine[2 + delta] * sine[2 + delta];
            const Rat near_factor = symmetric_sine * cosine[2 + delta];
            REQUIRE(near_factor > 0);
            const auto near_quadratic = scale(near, 1 / (drag_scale * near_factor));
            const Rat far_quadratic = far / (drag_scale * symmetric_sine);
            check_within(
                near_quadratic,
                rational_word(
                    condition.span_reference->at("quadratic_near_drag_N_per_rad2").get<double>()),
                4 * force_tolerance / near_factor);
            check_within(
                {far_quadratic, Rat{0}},
                rational_word(condition.span_reference->at("quadratic_trefftz_drag_N_per_rad2")
                                  .get<double>()),
                4 * force_tolerance / symmetric_sine);
        }
    }
}
} // namespace

TEST_CASE("A frozen finite wing reproduces its complete steady polar and loads",
          "[V-AERO][AE-02][UW-D03][UW-D05][N4-06][verification]") {
    std::ifstream source{UVLAS_FINITE_WING_CASE_FIXTURE};
    std::ifstream observations{UVLAS_FINITE_WING_STEADY_FIXTURE};
    REQUIRE(source);
    REQUIRE(observations);
    const auto frozen = nlohmann::json::parse(source);
    const auto expected = nlohmann::json::parse(observations);
    REQUIRE(frozen.at("cases").size() == 2);
    REQUIRE(expected.at("cases").size() == 2);
    constexpr std::size_t wing_index = UVLAS_FINITE_WING_INDEX;
    static_assert(wing_index < 2);
    const auto& wing = frozen.at("cases").at(wing_index);
    const auto& row = expected.at("cases").at(wing_index);
    CAPTURE(wing.at("dihedral_degrees").get<int>());
    REQUIRE(wing.at("dihedral_degrees") == row.at("dihedral_degrees"));

    constexpr int series = UVLAS_FINITE_WING_SERIES;
    static_assert(series >= 0 && series <= 2);
    if constexpr (series == 0) {
        std::ifstream polar_file{UVLAS_FINITE_WING_POLARS_FIXTURE};
        REQUIRE(polar_file);
        const auto polars = nlohmann::json::parse(polar_file);
        REQUIRE(polars.at("cases").size() == 2);
        REQUIRE(polars.at("band_bounds_m").size() == 21);
        for (std::size_t i = 0; i <= band_count; ++i) {
            const Rat boundary = band_border(i);
            CHECK(abs(rational_word(polars.at("band_bounds_m").at(i).get<double>()) - boundary) <
                  rational_word(1e-15));
        }
        const auto& polar = polars.at("cases").at(wing_index);
        REQUIRE(wing.at("dihedral_degrees") == polar.at("dihedral_degrees"));
        std::ifstream line_file{UVLAS_FINITE_WING_SPAN_LINES_FIXTURE};
        REQUIRE(line_file);
        const auto lines = nlohmann::json::parse(line_file);
        REQUIRE(lines.at("cases").size() == 2);
        const auto& line = lines.at("cases").at(wing_index);
        REQUIRE(wing.at("dihedral_degrees") == line.at("dihedral_degrees"));
        REQUIRE(line.at("angles").size() == 5);
        REQUIRE(line.at("slopes").size() == 2);
        check_line_invariants();
        run_wing(wing, row, expected, polar,
                 {frozen.at("wake_length_m").get<double>(),
                  frozen.at("core_radius_m").get<double>(), true, nullptr, &line});
    } else {
        std::ifstream sensitivity_file{UVLAS_FINITE_WING_SENSITIVITY_FIXTURE};
        REQUIRE(sensitivity_file);
        const auto sensitivity = nlohmann::json::parse(sensitivity_file);
        REQUIRE(sensitivity.at("cases").size() == 2);
        const auto& selected = sensitivity.at("cases").at(wing_index);
        REQUIRE(wing.at("dihedral_degrees") == selected.at("dihedral_degrees"));
        constexpr std::array closure_lengths{256., 512., 1024.};
        constexpr std::array core_radii{1e-4, 1e-5};
        const auto& ideal = selected.at("ideal_same_finite_closure");
        REQUIRE(ideal.at("wake_length_m") == 512.);
        REQUIRE(ideal.at("core_radius_m") == 0.);
        REQUIRE(ideal.at("angles").size() == 1);
        const auto& ideal_angle = ideal.at("angles").at(0);
        REQUIRE(ideal_angle.at("alpha_degrees") == 1);
        const auto& samples = selected.at(series == 1 ? "closure" : "core");
        REQUIRE(samples.size() == (series == 1 ? closure_lengths.size() : core_radii.size()));
        for (std::size_t i = 0; i < samples.size(); ++i) {
            const auto& sample = samples.at(i);
            const double length = series == 1 ? closure_lengths.at(i) : 512.;
            const double radius = series == 1 ? 1e-3 : core_radii.at(i);
            CAPTURE(length, radius);
            REQUIRE(sample.at("wake_length_m") == length);
            REQUIRE(sample.at("core_radius_m") == radius);
            run_wing(wing, sample, expected, sample,
                     {length, radius, false, length == 512. ? &ideal_angle : nullptr});
        }
    }
}
