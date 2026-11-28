#include "aero/closed_contour.hpp"
#include "aero/geometry_snapshot_fixture.hpp"
#include "aero/p1_shared_action.hpp"
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
#include <bit>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numbers>
#include <numeric>
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
    explicit FiniteWingMesh(std::span<const std::array<double, 3>> coordinates,
                            std::size_t chord_count = 2, std::size_t span_count = 4) {
        const auto point = [span_count](std::size_t chord, std::size_t span) {
            return a::PointId{1 + (span_count + 1) * chord + span};
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
        for (std::size_t chord = 0; chord <= chord_count; ++chord) {
            for (std::size_t span = 0; span < span_count; ++span) {
                const auto index = span_count * chord + span;
                const a::PhysicalFaceId face{100 + index};
                const a::BandId band{50 + span};
                const std::array vertex{point(chord, span), point(chord, span + 1),
                                        point(chord + 1, span + 1), point(chord + 1, span)};
                const std::array outline{edge(vertex[0], vertex[1]), edge(vertex[1], vertex[2]),
                                         edge(vertex[2], vertex[3]), edge(vertex[3], vertex[0])};
                faces.push_back({face, band});
                roles.push_back(
                    {face, chord == chord_count ? a::AeroFaceRole::Wake : a::AeroFaceRole::Bound});
                quads.emplace_back(
                    face, vertex, outline, edge(vertex[0], vertex[2]),
                    std::array{a::FacetId{200 + 2 * index}, a::FacetId{201 + 2 * index}});
                if (chord == chord_count)
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
    std::vector<Point> vertices;
    std::vector<Point> collocation;
    std::vector<cf::ExplicitCollocationTangents> tangents;
};

CaseLattice lattice(const io::SiCase& si, double wake_length, core::ResourcePlan& plan,
                    core::Index chord_count = 2, core::Index half_span = 2) {
    const auto root_case = take(take(si.source_case()).get().root());
    const auto model = root_case.model.record();
    const structure::PgReferenceLimits bounds{exact_limits(), 1'000'000, 1'000'000, {}};
    const auto make = [&](std::string_view number) {
        return take(Exact::from_decimal(number, exact_limits(), plan));
    };
    const auto center = make("0.5"), quarter = make("0.25"), three_quarters = make("0.75");
    const auto span_count = 2 * half_span;
    const auto vertex_width = span_count + 1;
    const auto panel_count = static_cast<std::size_t>(chord_count * span_count);
    CaseLattice result;
    result.vertices.resize(static_cast<std::size_t>((chord_count + 2) * vertex_width));
    result.collocation.resize(panel_count);
    result.tangents.resize(panel_count, {xyz({}), xyz({})});
    REQUIRE(model.aero_surfaces.size() == 2);
    REQUIRE(model.geometry.size() == 2);
    for (core::Index side = 0; side < 2; ++side) {
        const auto shape = take(model.aero_surfaces.at(side)).record();
        REQUIRE(std::holds_alternative<io::schema::AeroAttachmentFixed>(shape.attachment.record()));
        REQUIRE(shape.chord.size() == static_cast<std::uint64_t>(chord_count + 1));
        REQUIRE(shape.span.size() == static_cast<std::uint64_t>(half_span + 1));
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
        for (core::Index c = 0; c <= chord_count; ++c) {
            const auto first = chord(c);
            const auto xi = c == chord_count ? take(first.clone(exact_limits(), plan))
                                             : mix(first, chord(c + 1), quarter, plan);
            for (core::Index s = 0; s <= half_span; ++s) {
                const auto location = project(sample(xi, span(s)).position, plan);
                const auto index =
                    static_cast<std::size_t>(vertex_width * c + half_span * side + s);
                if (side == 1 && s == 0) {
                    CHECK(location == result.vertices[index]);
                } else {
                    result.vertices[index] = location;
                }
            }
        }
        for (core::Index c = 0; c < chord_count; ++c) {
            const auto xi = mix(chord(c), chord(c + 1), three_quarters, plan);
            for (core::Index s = 0; s < half_span; ++s) {
                const auto point = sample(xi, mix(span(s), span(s + 1), center, plan));
                const auto index = static_cast<std::size_t>(span_count * c + half_span * side + s);
                result.collocation[index] = project(point.position, plan);
                result.tangents[index] = {xyz(project(point.chord_tangent, plan)),
                                          xyz(project(point.span_tangent, plan))};
            }
        }
    }
    for (core::Index i = 0; i < vertex_width; ++i) {
        const auto tail = static_cast<std::size_t>((chord_count + 1) * vertex_width + i);
        result.vertices[tail] =
            result.vertices[static_cast<std::size_t>(chord_count * vertex_width + i)];
        result.vertices[tail][0] += wake_length;
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

template <class Add>
void visit_trefftz_midpoint_terms(const cf::DiscreteCirculation& solved, const FiniteWingMesh& mesh,
                                  const a::GlobalContactState& contact,
                                  std::span<const Point> vertices, std::size_t bound_count,
                                  std::size_t strips, Add add_term) {
    std::vector<RationalPoint> edge(strips + 1);
    std::vector<Rat> trace(strips);
    for (std::size_t i = 0; i < strips; ++i) {
        const auto wake =
            take_global(contact.uvlm_quad(mesh.quads[bound_count + i].facets[0], 1000000));
        const auto& owner = wake.attachment;
        if (!owner)
            FAIL("The Trefftz trace must come from a qualified trailing owner");
        const auto first = endpoint(vertices, wake.corners[0]);
        if (i != 0)
            CHECK(edge[i] == first);
        edge[i] = first;
        edge[i + 1] = endpoint(vertices, wake.corners[1]);
        const auto row = bound_count - strips + i;
        REQUIRE(take(solved.row(row)).face == owner->bound_face);
        trace[i] = exact_word(words(solved.solution(), row)) * owner->coefficient;
    }
    std::vector<Rat> jumps(strips + 1);
    for (std::size_t i = 0; i < strips; ++i) {
        jumps[i] -= trace[i];
        jumps[i + 1] += trace[i];
    }
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
            add_term(Rat{trace[i] * jumps[j] * (ry * dy + rz * dz) / r2});
        }
    }
}

class DyadicEnclosureSum {
public:
    void add(const Rat& term) {
        const Big scaled = numerator(term) * units_;
        const Big divisor = denominator(term);
        Big integral, remainder;
        boost::multiprecision::divide_qr(scaled, divisor, integral, remainder);
        integral_ += integral;
        // Integer division truncates toward zero. Each discarded remainder has
        // magnitude below one unit, irrespective of the sign of the term.
        if (remainder != 0)
            ++remainder_count_;
    }

    [[nodiscard]] Enclosure enclosure() const {
        return {Rat{integral_} / units_, Rat{remainder_count_} / units_};
    }

private:
    const Big units_ = Big{1} << 100;
    Big integral_{};
    std::size_t remainder_count_{};
};

Rat trefftz_midpoint(const cf::DiscreteCirculation& solved, const FiniteWingMesh& mesh,
                     const a::GlobalContactState& contact, std::span<const Point> vertices,
                     const Rat& density, std::size_t bound_count = columns,
                     std::size_t strips = 4) {
    Rat sum{};
    DyadicEnclosureSum enclosed_sum;
    visit_trefftz_midpoint_terms(solved, mesh, contact, vertices, bound_count, strips,
                                 [&](const Rat& term) {
                                     sum += term;
                                     enclosed_sum.add(term);
                                 });
    const auto enclosure = enclosed_sum.enclosure();
    CHECK(abs(sum - enclosure.center) <= enclosure.radius);
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

[[maybe_unused]] void run_wing(const nlohmann::json& wing, const nlohmann::json& expected,
                               const nlohmann::json& reference_metadata,
                               const nlohmann::json& polar, const WingCondition& condition) {
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
#ifdef UVLAS_FINITE_WING_REFINEMENT_FIXTURE
#ifdef UVLAS_FINITE_WING_FORMAL_PRODUCT
constexpr bool formal_product = true;
#else
constexpr bool formal_product = false;
#endif

#ifdef UVLAS_FINITE_WING_FORMAL_GEOMETRY
void verify_formal_wing_forecast(const CaseLattice& grid, const FiniteWingMesh& mesh,
                                 const a::OwnedGeometry& geometry,
                                 const a::GlobalContactState& contact,
                                 const a::GeometryLimits& geometry_caps, std::size_t span,
                                 a::GeometryWorkspace& workspace, core::ResourcePlan& memory) {
    const auto count = grid.collocation.size();
    std::vector<a::P1FaceStrengths> maps;
    std::vector<a::P1StrengthId> ids;
    maps.reserve(mesh.faces.size());
    ids.reserve(mesh.faces.size() + span);
    for (std::size_t i = 0; i < count; ++i) {
        const auto face = mesh.faces[i].id;
        maps.emplace_back(face, a::P0BoundColumn{});
        ids.push_back(a::P1StrengthId::bound(face));
    }
    for (std::size_t i = 0; i < span; ++i) {
        const auto panel = take_global(contact.uvlm_quad(mesh.quads[count + i].facets[0], 1000000));
        const auto& attachment = panel.attachment;
        if (!attachment)
            FAIL("Formal wing wake strip requires a trailing-edge attachment");
        const auto leading = attachment->coefficient == -1
                                 ? a::P1StrengthId{10000 + i}
                                 : a::P1StrengthId::bound(attachment->bound_face);
        const auto trailing = a::P1StrengthId{5000 + i};
        maps.emplace_back(panel.face, std::array{leading, trailing});
        ids.push_back(trailing);
        if (!leading.bound_face())
            ids.push_back(leading);
    }
    std::ranges::sort(ids);
    REQUIRE(maps.size() == (128 + 1) * 192);
    constexpr std::uint64_t query_cap = UINT64_C(150000000000);
    const a::P1AssemblyLimits admission{
        geometry_caps,
        {ls::EstimateAdmission{UINT64_C(8) * 1024 * 1024}, UINT64_C(1000000000)},
        query_cap};
    std::cout << "Formal geometry: points=" << mesh.ids.size() << " facets=" << 2 * maps.size()
              << " preparation_bytes=" << memory.allocated_bytes()
              << " kernel_step_cap=" << a::P1Workspace::maximum_steps() << '\n';
    for (const auto row : {std::size_t{0}, count / 2, count - 1}) {
        CAPTURE(row);
        const auto query = test_geometry::take(a::resolve_query(
            geometry, {a::QueryId{row + 1}, a::FreePoint{a::FreeSingle{grid.collocation[row]}}},
            geometry.context(), geometry_caps, workspace, memory));
        const a::P1ActionSource source{contact, query, {7, 11}};
        const auto forecast = test_geometry::take(
            a::forecast_bound_wake_serial_induction(source, maps, ids, admission));
        REQUIRE(forecast.reusable);
        CHECK(forecast.ordinary.reserved_steps <= query_cap);
        std::cout << "Formal forecast: row=" << row
                  << " preflight=" << forecast.ordinary.preflight_steps
                  << " reserved=" << forecast.ordinary.reserved_steps
                  << " first_bytes=" << forecast.initial_payload_bytes
                  << " repeated_bytes=" << forecast.repeated_payload_bytes
                  << " metadata_bytes=" << forecast.ordinary.estimated_metadata_bytes << '\n';
    }
}

void verify_formal_wing_contacts(const CaseLattice& grid, const FiniteWingMesh& mesh,
                                 const a::OwnedGeometry& geometry,
                                 const a::GlobalContactState& contact,
                                 const a::GeometryLimits& geometry_caps,
                                 a::GlobalLimits contact_caps, a::GeometryWorkspace& geometry_work,
                                 a::GlobalWorkspace& global_work, core::ResourcePlan& memory) {
    const auto index =
        take_global(a::ContourContactIndex::prepare(contact, contact_caps, global_work, memory));
    const auto query = [&](a::QueryDefinition definition) {
        return test_geometry::take(a::resolve_query(geometry, {a::QueryId{1}, definition},
                                                    geometry.context(), geometry_caps,
                                                    geometry_work, memory));
    };
    const auto qualify = [&](const a::ResolvedQuery& resolved, std::optional<a::EdgeId> self) {
        return index.qualify(resolved, self, geometry_caps, contact_caps, geometry_work,
                             global_work);
    };
    std::uint64_t collocation_candidates = 0, midpoint_candidates = 0;
    for (const auto& location : grid.collocation)
        collocation_candidates +=
            take_global(qualify(query(a::FreePoint{a::FreeSingle{location}}), {})).candidates;

    std::vector<a::EdgeId> midpoint_edges;
    for (const auto& quad : std::span{mesh.quads}.first(grid.collocation.size()))
        midpoint_edges.insert(midpoint_edges.end(), quad.edges.begin(), quad.edges.end());
    std::ranges::sort(midpoint_edges);
    const auto duplicates = std::ranges::unique(midpoint_edges);
    midpoint_edges.erase(duplicates.begin(), duplicates.end());
    for (const auto& trailing : mesh.trailing) {
        const auto found = std::ranges::lower_bound(midpoint_edges, trailing.edge);
        REQUIRE(found != midpoint_edges.end());
        REQUIRE(*found == trailing.edge);
        midpoint_edges.erase(found);
    }
    REQUIRE(midpoint_edges.size() == 49'280);
    for (const auto edge : midpoint_edges)
        midpoint_candidates +=
            take_global(qualify(query(a::SourceEdgePoint{edge, a::SingleWord{0.5}}), edge))
                .candidates;

    // The reference closed boxes used the exported Binary64 input vertices and
    // a one-ulp enclosing interval for each exact midpoint's coordinates.
    CHECK(collocation_candidates == 238'380);
    CHECK(midpoint_candidates >= 98'112);
    CHECK(midpoint_candidates <= 663'240);
    const auto on_vertex = qualify(query(a::FreePoint{a::FreeSingle{grid.vertices.front()}}), {});
    REQUIRE_FALSE(on_vertex);
    CHECK(on_vertex.error().origin.cause.code == core::FailureCode::GeometryInvalid);
    const auto diagonal = mesh.quads.front().diagonal;
    const auto wrong_self =
        qualify(query(a::SourceEdgePoint{diagonal, a::SingleWord{0.5}}), diagonal);
    REQUIRE_FALSE(wrong_self);
    CHECK(wrong_self.error().origin.cause.code == core::FailureCode::Unsupported);
    std::cout << "Formal contour contact: collocation_candidates=" << collocation_candidates
              << " midpoint_candidates=" << midpoint_candidates
              << " refit_steps=" << index.preparation_steps() << '\n';
}
#endif

void verify_prepared_contour(const FiniteWingMesh& mesh, const a::OwnedGeometry& geometry,
                             const a::GlobalContactState& contact,
                             std::span<const a::ResolvedQuery> queries,
                             std::span<const a::TrailingEdgeRow> transfers, a::EdgeId self,
                             const cf::DiscreteNonpenetrationLimits& caps,
                             a::GlobalLimits contact_caps, a::GeometryWorkspace& geometry_work,
                             a::GlobalWorkspace& global_work, a::SharedP1Workspace& shared,
                             a::LocalPacketWorkspace& local, a::P1Workspace& kernel,
                             core::ResourcePlan& memory) {
    std::vector<a::P1StrengthId> ids;
    std::vector<a::P1FaceStrengths> maps;
    for (const auto& face : std::span{mesh.faces}.first(queries.size())) {
        ids.push_back(a::P1StrengthId::bound(face.id));
        maps.emplace_back(face.id, a::P0BoundColumn{});
    }
    for (std::size_t i = 0; i < transfers.size(); ++i) {
        auto key = a::P1StrengthId::bound(transfers[i].bound_face);
        if (transfers[i].coefficient == -1) {
            key = a::P1StrengthId{i};
            ids.push_back(key);
        }
        maps.emplace_back(mesh.faces[queries.size() + i].id, a::ConstantWakeColumn{key});
    }
    std::ranges::sort(ids);
    const auto check_bounds = [](a::ContourInterval candidate, const a::ScaledEnclosure& exact) {
        const auto value = unpack(exact);
        CHECK(rational_word(candidate.lower) <= value.center - value.radius);
        CHECK(rational_word(candidate.upper) >= value.center + value.radius);
        CHECK(candidate.upper - candidate.lower < 1e-10);
    };
    core::ResourcePlan measured{memory.remaining_bytes()};
    const auto contour = test_geometry::take(
        a::ClosedContour::prepare(contact, maps, ids, contact_caps, global_work, measured));
    CHECK(contour.preparation_steps() > 0);
    core::ResourcePlan one_short{measured.allocated_bytes() - 1};
    const auto denied =
        a::ClosedContour::prepare(contact, maps, ids, contact_caps, global_work, one_short);
    REQUIRE_FALSE(denied);
    CHECK(denied.error().cause.code == core::FailureCode::ResourceExceeded);
    CHECK(one_short.allocated_bytes() == 0);
    auto workspace = test_geometry::take(contour.workspace(memory));
    std::vector<a::ContourInterval> row(ids.size());
    const a::ContourVector normal{{{}, {}, {1, 1}}};
    const auto& query = queries.front();
    const auto steps =
        test_geometry::take(contour.project(query, normal, row, caps.induction.geometry,
                                            contact_caps, geometry_work, global_work, workspace));
    auto too_short = contact_caps;
    too_short.max_steps = steps - 1;
    const auto exceeded = contour.project(query, normal, row, caps.induction.geometry, too_short,
                                          geometry_work, global_work, workspace);
    REQUIRE_FALSE(exceeded);
    CHECK(exceeded.error().cause.code == core::FailureCode::ResourceExceeded);
    CHECK(test_geometry::take(contour.project(query, normal, row, caps.induction.geometry,
                                              contact_caps, geometry_work, global_work,
                                              workspace)) == steps);

    constexpr a::P1ActionStamp action_stamp{125, 2};
    const a::P1ActionSource source{contact, query, action_stamp};
    const auto full = test_geometry::take(a::prepare_bound_wake_induction(
        source, maps, ids, caps.geometry_accuracy, caps.induction_accuracy, caps.induction, shared,
        local, kernel, memory));
    const a::P1ComplexWords zero{a::SingleWord{0}, a::SingleWord{0}};
    const a::P1ComplexWords one{a::SingleWord{1}, a::SingleWord{0}};
    std::vector<a::P1ComplexWords> strengths(ids.size(), zero);
    for (std::size_t i = 0; i < strengths.size(); ++i) {
        strengths[i] = one;
        const auto reference = test_geometry::take(
            a::apply_shared_p1(full, source, strengths, {}, caps.action_steps, shared));
        check_bounds(row[i], reference.velocity()[2].real);
        strengths[i] = zero;
    }

    const auto midpoint = test_geometry::take(
        a::resolve_query(geometry, {a::QueryId{812}, a::SourceEdgePoint{self, a::SingleWord{0.5}}},
                         geometry.context(), caps.induction.geometry, geometry_work, memory));
    const a::P1ActionSource midpoint_source{contact, midpoint, action_stamp};
    const auto midpoint_action = test_geometry::take(a::prepare_bound_wake_midpoint_induction(
        midpoint_source, self, maps, ids, caps.geometry_accuracy, caps.induction_accuracy,
        caps.induction, shared, local, kernel, memory));
    std::vector<a::ContourInterval> weights(ids.size());
    for (std::size_t i = 0; i < strengths.size(); ++i) {
        const double value = static_cast<double>(i) / 8. - .25;
        weights[i] = {value, value};
        strengths[i].real = a::SingleWord{value};
    }
    const auto candidate =
        test_geometry::take(contour.apply(midpoint, self, weights, caps.induction.geometry,
                                          contact_caps, geometry_work, global_work, workspace));
    const auto reference = test_geometry::take(a::apply_shared_p1(
        midpoint_action, midpoint_source, strengths, {}, caps.action_steps, shared));
    for (std::size_t axis = 0; axis < candidate.first.size(); ++axis)
        check_bounds(candidate.first[axis], reference.velocity()[axis].real);
}

#ifdef UVLAS_FINITE_WING_FORMAL_PRODUCT
void check_formal_product(const cf::DiscreteCirculationLoads& loads, const CaseLattice& grid,
                          const FiniteWingMesh& mesh, const a::OwnedGeometry& geometry,
                          const a::GlobalContactState& contact,
                          std::span<const a::TrailingEdgeRow> transfers, const Rat& sine,
                          const nlohmann::json& common_limit) {
    const auto count = grid.collocation.size();
    const auto ny = transfers.size();
    const Rat density = rational_word(1.2);
    const Rat sine_squared = sine * sine;
    const auto force = unpack_vector(loads.total_force());
    const auto axial = scale(force[0], Rat{1} / sine_squared);
    DyadicEnclosureSum trefftz_terms;
    visit_trefftz_midpoint_terms(loads.primal(), mesh, contact, grid.vertices, count, ny,
                                 [&](const Rat& term) { trefftz_terms.add(term); });
    const auto trefftz = scale(trefftz_terms.enclosure(),
                               -density / (4 * rational_word(std::numbers::pi) * sine_squared));
    // The streamwise ambient term is the odd part of Fz under incidence reversal.
    // Every closed bound ring has zero directed span displacement. With its trailing
    // edge omitted from the loads, the remaining displacement is minus that edge's.
    Rat trailing = 0;
    for (std::size_t k = 0; k < ny; ++k) {
        const auto panel = count - ny + k;
        const auto te = transfers[k].edge;
        const auto gamma = exact_word(words(loads.primal().solution(), panel));
        std::size_t sides = 0;
        for (const auto facet : mesh.quads[panel].facets) {
            const auto edges = test_geometry::take(geometry.facet_edges(facet));
            const auto turns = test_geometry::take(geometry.incidence(facet));
            for (std::size_t side = 0; side < edges.size(); ++side) {
                const auto edge = take_global(contact.uvlm_load_edge(edges[side], 1'000'000));
                if (edge.edge != te)
                    continue;
                const auto y0 = grid.vertices.at(edge.endpoints[0].value() - 1)[1];
                const auto y1 = grid.vertices.at(edge.endpoints[1].value() - 1)[1];
                trailing +=
                    gamma * static_cast<int>(turns[side]) * (rational_word(y1) - rational_word(y0));
                ++sides;
            }
        }
        REQUIRE(sides == 1);
    }
    const Rat near = axial.center + trailing * density * 20 / sine;
    const Rat limit = rational_word(common_limit.at("reference_limit").get<double>());
    const Rat uncertainty =
        rational_word(common_limit.at("reference_relative_uncertainty").get<double>());
    const Rat trefftz_lower = trefftz.center - trefftz.radius;
    REQUIRE(trefftz_lower > 0);
    REQUIRE(limit > 0);
    REQUIRE(uncertainty >= 0);
    const Rat gap = (abs(near - trefftz.center) + axial.radius + trefftz.radius) / trefftz_lower;
    const Rat indicator = (abs(near - limit) + axial.radius) / limit + uncertainty;
    const Rat threshold = Rat{1} / 100;
    CHECK(gap <= threshold);
    CHECK(indicator <= threshold);
    const Rat finite_trefftz = rational_word(common_limit.at("candidate_trefftz").get<double>());
    CHECK(abs(trefftz.center - finite_trefftz) + trefftz.radius < Rat{1} / 10'000'000);

    // Verify the assumed lower bound on |n_z| with the actual input tangents.
    const Rat minimum_normal = Rat{99} / 100;
    const Rat normal_squared = minimum_normal * minimum_normal;
    Rat residual_squared = 0;
    for (std::size_t i = 0; i < count; ++i) {
        const auto& tangents = grid.tangents[i];
        CHECK(exact_word(tangents[0][1]) == 0);
        CHECK(exact_word(tangents[0][2]) == 0);
        const auto transverse = exact_word(tangents[1][1]);
        const auto height = exact_word(tangents[1][2]);
        REQUIRE(transverse != 0);
        CHECK(normal_squared * (transverse * transverse + height * height) <=
              transverse * transverse);
        const auto row = unpack(take(loads.primal().row(i)).original_residual);
        const Rat bound = abs(row.center) + row.radius;
        residual_squared += bound * bound;
    }
    const Rat rhs_squared_lower = 400 * count * normal_squared * sine_squared;
    const Rat residual_ratio_squared = residual_squared / rhs_squared_lower;
    const Rat residual_budget = Rat{1} / UINT64_C(100000000000);
    CHECK(residual_ratio_squared <= residual_budget * residual_budget);
    std::cout << std::setprecision(17)
              << "Formal midpoint drag N/rad^2, near=" << near.convert_to<double>()
              << ", radius=" << axial.radius.convert_to<double>()
              << ", Trefftz=" << trefftz.center.convert_to<double>()
              << ", Trefftz radius=" << trefftz.radius.convert_to<double>()
              << ", near/far gap=" << gap.convert_to<double>()
              << ", ideal reference indicator=" << indicator.convert_to<double>()
              << ", relative L2 residual ceiling="
              << std::sqrt(residual_ratio_squared.convert_to<double>()) << '\n';
}
#endif

void verify_harmonic_grid(const CaseLattice& grid, const FiniteWingMesh& mesh,
                          a::GeometryInput geometry_input, const a::GeometryLimits& geometry_caps,
                          const a::GlobalLimits& contact_caps,
                          cf::DiscreteNonpenetrationLimits action_limits,
                          cf::DiscreteCirculationRequest request) {
    std::ifstream fixture(UVLAS_HARMONIC_M2N4_FIXTURE);
    REQUIRE(fixture);
    const auto reference = nlohmann::json::parse(fixture);
    REQUIRE(reference.at("schema") == "finite-wing-harmonic-m2n4/1");
    REQUIRE(reference.at("cases").size() == 6);
    const auto wake_length = reference.at("wake_length_m").get<double>();
    auto high = mesh.xyz;
    std::vector<double> low(high.size());
    constexpr std::size_t count = 8;
    constexpr std::size_t span = 4;
    constexpr std::size_t width = span + 1;
    REQUIRE(grid.vertices.size() == 4 * width);
    for (std::size_t i = 0; i < width; ++i) {
        const auto x = high[3 * (2 * width + i)];
        const auto tail = 3 * (3 * width + i);
        REQUIRE(std::abs(x) < wake_length);
        // FastTwoSum with the wake length as the larger operand. All four
        // bands recover the same exact translation from the retained words.
        high[tail] = wake_length + x;
        low[tail] = x - (high[tail] - wake_length);
        CHECK(rational_word(high[tail]) + rational_word(low[tail]) - rational_word(x) ==
              rational_word(wake_length));
    }
    const a::DoublePositions positions{high, low};
    geometry_input.reference = positions;
    geometry_input.current = positions;
    geometry_input.core_radius = a::SingleWord{reference.at("core_m").get<double>()};
    core::ResourcePlan plan{UINT64_C(2000000000)};
    auto workspace = test_geometry::workspace(plan);
    const auto geometry = test_geometry::take(
        a::prepare_owned_geometry(geometry_input, geometry_caps, *workspace, plan));
    auto global = take(a::GlobalWorkspace::create(plan, contact_caps.metadata));
    const auto topology = take_global(a::build_global_topology(
        geometry, {mesh.roles, mesh.owners, mesh.trailing}, contact_caps, *global, plan));
    const auto contact =
        take_global(a::check_global_geometry(topology, geometry, contact_caps, *global, plan));
    std::vector<a::ResolvedQuery> queries;
    queries.reserve(count);
    std::vector<cf::DiscreteCollocation> rows;
    rows.reserve(count);
    for (std::size_t i = 0; i < count; ++i) {
        queries.push_back(test_geometry::take(a::resolve_query(
            geometry, {a::QueryId{i + 1}, a::FreePoint{a::FreeSingle{grid.collocation[i]}}},
            geometry.context(), geometry_caps, *workspace, plan)));
        rows.emplace_back(mesh.quads[i].facets[0], queries.back(), grid.tangents[i], xyz({}),
                          a::SingleWord{0}, a::SingleWord{0});
    }
    std::vector<cf::DiscreteWakePanel> wake;
    for (std::size_t i = 0; i < span; ++i)
        wake.emplace_back(mesh.quads[count + i].facets[0], cf::AttachedWakeFront{},
                          cf::SteadyWakeClosure{});
    constexpr a::P1ActionStamp trial{815, 6};
    const cf::DiscreteNonpenetrationInput input{contact,
                                                rows,
                                                wake,
                                                {core::AcceptedStateId{2}, {}},
                                                trial,
                                                stamp,
                                                p::RequestedSide::Automatic};
    request.scales.context = context(geometry.context());
    request.certified_closed_contour = true;
    request.solver_limits.maximum_update_rank = span;
    request.matrix_absolute = 1e-10;
    request.algebraic.norm = ls::Norm::L2;
    action_limits.action_steps = action_limits.induction.steps;
    const cf::HarmonicLoadInput loading{a::SingleWord{reference.at("density_kg_m3").get<double>()},
                                        xyz({.375, 0., 0.})};
    const auto speed = reference.at("inflow_m_s").get<double>();
    std::vector<cf::HarmonicReflectionPair> reflected;
    reflected.reserve(count / 2);
    for (std::size_t i = 0; i < count / 2; ++i) {
        const auto first = i / (span / 2) * span + i % (span / 2);
        reflected.push_back(
            {mesh.faces[first].id, mesh.faces[first / span * span + span - 1 - first % span].id});
    }
    for (const auto& sample : reference.at("cases")) {
        const auto motion = sample.at("motion").get<std::string>();
        CAPTURE(sample.at("k").get<double>(), motion);
        REQUIRE((motion == "heave" || motion == "pitch"));
        const bool heave = motion == "heave";
        const auto starts = sample.at("panel_starts").get<std::vector<double>>();
        const auto omega = sample.at("omega").get<double>();
        const cf::StraightWakeHarmonicInput phase{xyz({speed, 0., 0.}), starts, omega,
                                                  16 * starts.size()};
        std::vector<cf::HarmonicCirculationRhs> forcing;
        forcing.reserve(count);
        double rhs_square = 0.;
        for (std::size_t i = 0; i < count; ++i) {
            const auto rhs = sample.at("rhs").at(i).get<std::array<double, 2>>();
            // The horizontal wing's normal component is one. The pivot is
            // the same x=3/8 m used for the output moment.
            const double real = heave ? 0. : -speed;
            const double imag = heave ? omega : -omega * (grid.collocation[i][0] - .375);
            CHECK(std::hypot(real - rhs[0], imag - rhs[1]) <= 1e-13 * std::hypot(rhs[0], rhs[1]));
            forcing.push_back({mesh.faces[i].id, a::SingleWord{real}, a::SingleWord{imag}});
            rhs_square += real * real + imag * imag;
        }
        request.original_absolute = .99 * std::sqrt(rhs_square) * 1e-11 / 20.;
        for (const bool condensed : {false, true}) {
            CAPTURE(condensed);
            request.contour_workers = condensed ? 2 : 1;
            auto selected = phase;
            if (condensed)
                selected.reflection_pairs = reflected;
            const auto result = take(cf::solve_discrete_harmonic_loads(
                input, selected, forcing, loading, action_limits, request, plan));
            REQUIRE(result.loads);
            CHECK(result.expanded.has_value() == condensed);
            CHECK(result.physical_values().shape().rows == static_cast<core::Index>(count));
            CHECK(result.solution.values().shape().rows ==
                  static_cast<core::Index>(condensed ? count / 2 : count));
            CHECK(result.residual.size() == count);
            const auto& solution = std::get<numerics::TwoWord<std::complex<double>>>(
                result.physical_values().storage());
            for (std::size_t i = 0; i < count; ++i) {
                const auto expected = sample.at("gamma").at(i).get<std::array<double, 2>>();
                const std::complex<double> gamma{expected[0], expected[1]};
                CHECK(std::abs(solution.hi.values()[i] + solution.lo.values()[i] - gamma) <=
                      1e-9 * std::abs(gamma));
            }
            const auto check = [&](std::size_t output, const a::ScaledEnclosure& real,
                                   const a::ScaledEnclosure& imag) {
                const auto expected = sample.at("loads").at(output).get<std::array<double, 2>>();
                const auto tolerance = rational_word(5e-10 * std::hypot(expected[0], expected[1]));
                check_within(unpack(real), rational_word(expected[0]), tolerance);
                check_within(unpack(imag), rational_word(expected[1]), tolerance);
            };
            check(0, result.loads->force_real[2], result.loads->force_imaginary[2]);
            check(1, result.loads->moment_real[1], result.loads->moment_imaginary[1]);
        }
        if (heave && sample.at("k").get<double>() == .05) {
            auto selected = phase;
            selected.reflection_pairs = reflected;
            forcing[0].real = a::SingleWord{1.};
            const auto asymmetric = cf::solve_discrete_harmonic_circulation(
                input, selected, forcing, action_limits, request, plan);
            REQUIRE_FALSE(asymmetric);
            CHECK(asymmetric.error().code == core::FailureCode::InvalidInput);
            forcing[0].real = a::SingleWord{0.};
            auto wrong = reflected;
            std::swap(wrong[0][1], wrong[1][1]);
            selected.reflection_pairs = wrong;
            const auto wrong_geometry = cf::solve_discrete_harmonic_circulation(
                input, selected, forcing, action_limits, request, plan);
            REQUIRE_FALSE(wrong_geometry);
            CHECK(wrong_geometry.error().code == core::FailureCode::InvalidInput);
        }
    }
}

void verify_harmonic_thin_panel() {
    struct Example {
        std::array<Point, 4> vertices;
        Point query;
        double circulation;
    };
    // Swept panels and exact Binary64 collocations from formal M256/N192 rows
    // (0, 162) and (255, 175), using zero frequency to isolate the closed contour.
    // Reference circulation is the reciprocal of the directed six-edge contour
    // influence, evaluated from exact binary inputs with a 100-digit antiderivative.
    const std::array examples{
        Example{{Point{2.1784515365038395, 3.1111404660392044, 0.},
                 Point{2.216133040577926, 3.164955393735604, 0.},
                 Point{2.178496713893399, 3.1111404660392044, 0.},
                 Point{2.2161775341817913, 3.164955393735604, 0.}},
                {2.1973072340815776, 3.138047929887404, 0.},
                -7.740080467966106e-5},
        Example{{Point{3.2574907744621195, 3.698404363053158, 0.},
                 Point{3.2734800322953057, 3.732050807568877, 0.},
                 Point{3.2575096326705406, 3.698404363053158, 0.},
                 Point{3.2734986767383654, 3.732050807568877, 0.}},
                {3.26549790426254, 3.7152275853110175, 0.},
                -7.145052610271387e-5},
    };
    for (const auto& example : examples) {
        CAPTURE(example.query);
        auto tail_left = example.vertices[2];
        auto tail_right = example.vertices[3];
        tail_left[0] += 128.;
        tail_right[0] += 128.;
        const std::array vertices{example.vertices[0], example.vertices[1], example.vertices[2],
                                  example.vertices[3], tail_left,           tail_right};
        const FiniteWingMesh mesh{vertices, 1, 1};
        std::array<double, 18> low{};
        low[12] = vertices[2][0] - (tail_left[0] - 128.);
        low[15] = vertices[3][0] - (tail_right[0] - 128.);
        const a::DoublePositions positions{mesh.xyz, low};
        const std::array bands{a::BandId{50}};
        const auto geometry_caps = test_geometry::limits();
        const a::GlobalLimits contact_caps{std::numeric_limits<std::uint64_t>::max() / 4,
                                           ls::EstimateAdmission{65536}};
        core::ResourcePlan plan{UINT64_C(2000000000)};
        auto workspace = test_geometry::workspace(plan);
        const Number instant{5, 4};
        const auto geometry =
            test_geometry::take(a::prepare_owned_geometry({mesh.ids,
                                                           mesh.bindings,
                                                           positions,
                                                           positions,
                                                           mesh.edges,
                                                           mesh.faces,
                                                           bands,
                                                           {},
                                                           mesh.quads,
                                                           a::SingleWord{1e-6},
                                                           test_geometry::context(),
                                                           instant.view()},
                                                          geometry_caps, *workspace, plan));
        auto global = take(a::GlobalWorkspace::create(plan, contact_caps.metadata));
        const auto topology = take_global(a::build_global_topology(
            geometry, {mesh.roles, mesh.owners, mesh.trailing}, contact_caps, *global, plan));
        const auto contact =
            take_global(a::check_global_geometry(topology, geometry, contact_caps, *global, plan));
        const auto query = test_geometry::take(
            a::resolve_query(geometry, {a::QueryId{1}, a::FreePoint{a::FreeSingle{example.query}}},
                             geometry.context(), geometry_caps, *workspace, plan));
        const std::array rows{cf::DiscreteCollocation{
            mesh.quads[0].facets[0], query,
            cf::ExplicitCollocationTangents{xyz({1., 0., 0.}), xyz({0., 1., 0.})}, xyz({}),
            a::SingleWord{0.}, a::SingleWord{0.}}};
        const std::array wake{cf::DiscreteWakePanel{
            mesh.quads[1].facets[0], cf::AttachedWakeFront{}, cf::SteadyWakeClosure{}}};
        const cf::DiscreteNonpenetrationInput input{contact,
                                                    rows,
                                                    wake,
                                                    {core::AcceptedStateId{2}, {}},
                                                    a::P1ActionStamp{815, 6},
                                                    stamp,
                                                    p::RequestedSide::Automatic};
        const std::array rhs{
            cf::HarmonicCirculationRhs{mesh.faces[0].id, a::SingleWord{1.}, a::SingleWord{0.}}};
        const std::array starts{0.};
        const cf::StraightWakeHarmonicInput phase{xyz({20., 0., 0.}), starts, 0., 16};
        const std::array scales{20.};
        const std::array<ls::UnitDimensions, 1> gamma_units{{{2, 0, -1, 0, 0, 0, 0}}};
        const std::array<ls::UnitDimensions, 1> velocity_units{{{1, 0, -1, 0, 0, 0, 0}}};
        cf::DiscreteCirculationRequest request{
            {scales, scales, gamma_units, velocity_units, context(geometry.context())},
            {1, 1, 1, 5, ls::EstimateAdmission{1U << 20U}},
            {ls::Norm::L2, 1e-12, 0., 1e-12, core::ObservationId{20}, core::ObservationId{21}},
            core::RhsGeneration{3},
            core::NumericalGeneration{4},
            6e-9,
            1e-13,
            4.9e-13,
            0.};
        request.certified_closed_contour = true;
        const cf::DiscreteNonpenetrationLimits action_limits{
            {core::Precision::TwoWord, 1e-12},
            {core::Precision::TwoWord, 1e-12, 1e-12},
            {geometry_caps, {ls::EstimateAdmission{65536}, 1000000000}, UINT64_C(150000000000)},
            UINT64_C(150000000000),
            UINT64_C(70000000000000)};
        const auto solve = [&] {
            return cf::solve_discrete_harmonic_circulation(input, phase, rhs, action_limits,
                                                           request, plan);
        };
        auto initial = solve();
        INFO("Operation " << (initial ? -1 : static_cast<int>(initial.error().operation)));
        const auto qualified = take(std::move(initial));
        REQUIRE(qualified.residual.size() == 1);
        REQUIRE(qualified.residual.values().front());
        const auto& words = std::get<numerics::TwoWord<std::complex<double>>>(
            qualified.physical_values().storage());
        const auto gamma = words.hi.values().front() + words.lo.values().front();
        CHECK(std::abs(gamma - example.circulation) <= 2e-12 * std::abs(example.circulation));
        const auto original = unpack(qualified.residual.values().front()->real);
        CHECK(abs(original.center) + original.radius < rational_word(1e-11));
        request.matrix_absolute = 1e-16;
        const auto too_tight = solve();
        REQUIRE_FALSE(too_tight);
        CHECK(too_tight.error().code == core::FailureCode::PrecisionUnresolved);
    }
}

void run_refined_drag(const nlohmann::json& wing, const nlohmann::json& polar,
                      const nlohmann::json& oracle, double length, double radius) {
    const auto m = oracle.at("chord_panels").get<core::Index>();
    const auto n = oracle.at("span_panels").get<core::Index>();
    CAPTURE(m, n);
#ifdef UVLAS_FINITE_WING_FORMAL_GEOMETRY
    REQUIRE((m == 128 && n == 192));
#else
    REQUIRE(((m == 2 && n == 4) || (m == 4 && n == 8) || (m == 4 && n == 12)));
#endif
    const auto nx = static_cast<std::size_t>(m);
    const auto ny = static_cast<std::size_t>(n);
    const auto order = m * n;
    const auto count = static_cast<std::size_t>(order);
    auto recipe = finite_wing_case(wing, radius);
    auto& surfaces = recipe.field("model").field("aero_surfaces");
    for (auto& surface : surfaces.elements) {
        surface.set("chord", parse(oracle.at("chord_fractions").dump()));
        surface.set("span", parse(oracle.at("half_span_fractions").dump()));
    }
    const auto si = materialize(recipe, io::SourceId{534});
    core::ResourcePlan memory{m == 128 ? UINT64_C(24) * (UINT64_C(1) << 30) : UINT64_C(2000000000)};
    const auto grid = lattice(si, length, memory, m, n / 2);
    for (std::size_t i = 0; i < grid.vertices.size(); ++i) {
        const auto expected = oracle.at("vertices_m").at(i).get<Point>();
        for (std::size_t axis = 0; axis < expected.size(); ++axis)
            CHECK(std::abs(grid.vertices[i][axis] - expected[axis]) <= 2e-14);
    }
    for (std::size_t i = 0; i < count; ++i) {
        const auto expected = oracle.at("collocations_m").at(i).get<Point>();
        for (std::size_t axis = 0; axis < expected.size(); ++axis)
            CHECK(std::abs(grid.collocation[i][axis] - expected[axis]) <= 2e-14);
    }
    const FiniteWingMesh mesh{grid.vertices, nx, ny};
    std::vector<a::BandId> bands;
    bands.reserve(ny);
    for (std::size_t i = 0; i < ny; ++i)
        bands.emplace_back(50 + i);
    auto geometry_caps = test_geometry::limits();
    geometry_caps.counts = {mesh.ids.size(), mesh.edges.size(),     mesh.faces.size(),
                            bands.size(),    2 * mesh.faces.size(), mesh.quads.size()};
    geometry_caps.max_steps = 20000000;
    geometry_caps.predicate.max_steps = 1000000;
    auto geometry_work = test_geometry::workspace(memory);
    const Number instant{5, 4};
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
    REQUIRE(topology.trailing_edges().size() == ny);

#ifdef UVLAS_FINITE_WING_FORMAL_GEOMETRY
    if (m == 128 && !formal_product) {
        verify_formal_wing_forecast(grid, mesh, geometry, contact, geometry_caps, ny,
                                    *geometry_work, memory);
        verify_formal_wing_contacts(grid, mesh, geometry, contact, geometry_caps, contact_caps,
                                    *geometry_work, *global, memory);
        return;
    }
#endif

    const cf::DiscreteNonpenetrationLimits caps{
        {core::Precision::TwoWord, 1e-12},
        {core::Precision::TwoWord, 1e-12, 1e-12},
        {geometry_caps, {ls::EstimateAdmission{65536}, 1000000000}, UINT64_C(150000000000)},
        20000000,
        UINT64_C(70000000000000)};
    std::vector<a::ResolvedQuery> queries;
    queries.reserve(count);
    for (std::size_t i = 0; i < count; ++i)
        queries.push_back(test_geometry::take(a::resolve_query(
            geometry, {a::QueryId{i + 1}, a::FreePoint{a::FreeSingle{grid.collocation[i]}}},
            geometry.context(), geometry_caps, *geometry_work, memory)));
    std::vector<cf::DiscreteCollocation> rows;
    rows.reserve(count);
    std::vector<cf::BoundCirculationRate> rates;
    rates.reserve(count);
    std::vector<a::EdgeId> physical_edges;
    physical_edges.reserve(4 * count);
    for (std::size_t i = 0; i < count; ++i) {
        rows.emplace_back(mesh.quads[i].facets[0], queries[i], grid.tangents[i], xyz({}),
                          a::SingleWord{0}, a::SingleWord{0});
        rates.emplace_back(mesh.faces[i].id, a::SingleWord{0});
        const auto& edges = mesh.quads[i].edges;
        physical_edges.insert(physical_edges.end(), edges.begin(), edges.end());
    }
    std::ranges::sort(physical_edges);
    const auto duplicate = std::ranges::unique(physical_edges);
    physical_edges.erase(duplicate.begin(), duplicate.end());
    std::vector<cf::EdgeMidpointVelocity> midpoint;
    midpoint.reserve(nx * (2 * ny + 1));
    for (const auto edge : physical_edges)
        if (!std::ranges::any_of(mesh.trailing,
                                 [edge](const auto& row) { return row.edge == edge; }))
            midpoint.emplace_back(edge, xyz({}));
    REQUIRE(midpoint.size() == nx * (2 * ny + 1));
    std::vector<cf::DiscreteWakePanel> wake, frozen_wake;
    std::vector<a::TrailingEdgeRow> transfers;
    wake.reserve(ny);
    frozen_wake.reserve(ny);
    transfers.reserve(ny);
    for (std::size_t i = 0; i < ny; ++i) {
        const auto facet = mesh.quads[count + i].facets[0];
        wake.emplace_back(facet, cf::AttachedWakeFront{}, cf::SteadyWakeClosure{});
        frozen_wake.emplace_back(facet, cf::AttachedWakeFront{}, a::P1StrengthId{700 + i});
        const auto panel = take_global(contact.uvlm_quad(facet, 1'000'000));
        const auto& attachment = panel.attachment;
        if (!attachment)
            FAIL("The fixed steady strip requires a trailing-edge transfer");
        CHECK(attachment->bound_face == mesh.faces[count - ny + i].id);
        transfers.push_back(*attachment);
    }

    const std::vector<double> scales(count, 20.);
    const std::vector<ls::UnitDimensions> gamma_units(count, {2, 0, -1, 0, 0, 0, 0});
    const std::vector<ls::UnitDimensions> speed_units(count, {1, 0, -1, 0, 0, 0, 0});
    const cf::DiscreteCirculationRequest request{
        {scales, scales, gamma_units, speed_units, context(geometry.context())},
        {order, 1, 0, 5,
         ls::EstimateAdmission{formal_product ? UINT64_C(14) * (UINT64_C(1) << 30U) : 1U << 20U}},
        {formal_product ? ls::Norm::L2 : ls::Norm::Infinity, formal_product ? 1e-12 : 1e-25, 0.,
         formal_product ? 1e-12 : 1e-23, core::ObservationId{20}, core::ObservationId{21}},
        core::RhsGeneration{3},
        core::NumericalGeneration{4},
        1e-13,
        1e-13,
        1e-11,
        0.};
    constexpr a::P1ActionStamp trial{815, 6};
    if (count == 8 && wing.at("dihedral_degrees").get<int>() == 0) {
        verify_harmonic_grid(grid, mesh, input_geometry, geometry_caps, contact_caps, caps,
                             request);
        verify_harmonic_thin_panel();
    }
    const cf::DiscreteVectorLoadInput write_loads{
        trial, request.retained_generation, midpoint,
        rates, a::SingleWord{1.2},          xyz({.375, 0., 0.})};
    a::UniformFlowWorkspace flow_work;
    a::SharedP1Workspace shared;
    auto local = take(a::LocalPacketWorkspace::create(memory));
    auto kernel = take(a::P1Workspace::create(memory));
    if (count == 8) {
        core::ResourcePlan contour_test{UINT64_C(2000000000)};
        verify_prepared_contour(mesh, geometry, contact, queries, transfers, midpoint.front().edge,
                                caps, contact_caps, *geometry_work, *global, shared, *local,
                                *kernel, contour_test);
    }
    const auto gamma_prime = oracle.find("gamma_derivative_m2_s_per_rad");
    if (!formal_product) {
        REQUIRE(gamma_prime != oracle.end());
        REQUIRE(gamma_prime->size() == count);
    }
    std::array<Enclosure, 2> near{};
    std::array<Rat, 2> far{}, sin{}, cos{};
    for (std::size_t sample = 0; sample < 2; ++sample) {
        const auto& angle = polar.at("angles").at(1 + 2 * sample);
        REQUIRE(angle.at("alpha_degrees") == (sample == 0 ? -1 : 1));
        auto angled = recipe;
        angled.field("flow").field("uniform").set("value", parse(angle.at("flow_m_s").dump()));
        const auto angled_si = materialize(angled, io::SourceId{534});
        const auto field =
            take(cf::prepare_case_uniform_flow(angled_si, stamp, {}, flow_work, limits, memory));
        const auto project_flow = [&] {
            auto evaluated = take(cf::evaluate_case_uniform_flow(
                field, field, {stamp, instant.view()}, flow_work, limits.flow, memory));
            const Number one{1}, absolute{1, 1000000000}, zero{0};
            const cf::CaseFlowComponentBudget component{one.view(), absolute.view(), zero.view()};
            const std::array bounds{component, component, component};
            return take(cf::project_case_uniform_flow(
                std::move(evaluated), {core::Precision::TwoWord, bounds, bounds, {4096, budget}}));
        };
        auto flow = project_flow();
        sin[sample] = rational_from_view(take(flow.source().jet.velocity(2))) / 20;
        cos[sample] = rational_from_view(take(flow.source().jet.velocity(0))) / 20;

        const cf::DiscreteNonpenetrationInput input{contact,
                                                    rows,
                                                    wake,
                                                    {core::AcceptedStateId{2}, {}},
                                                    trial,
                                                    stamp,
                                                    p::RequestedSide::Automatic};
        auto selected = request;
        selected.certified_closed_contour = count == 8 && sample == 0;
        if (formal_product) {
            REQUIRE(rational_from_view(take(flow.source().jet.velocity(1))) == 0);
            selected.certified_closed_contour = true;
            selected.contour_workers = 4;
            selected.solver_limits.dense_residual_workers = 4;
            selected.matrix_absolute = 6e-9;
            selected.original_absolute = .99 * std::abs(sin[sample].convert_to<double>()) *
                                         std::sqrt(static_cast<double>(count)) * 1e-11;
        }
        auto selected_caps = caps;
        // The contour action also pays for exact contact qualification. Preserve
        // the fixture's existing induction ceiling for the positive contour run.
        if (selected.certified_closed_contour)
            selected_caps.action_steps = caps.induction.steps;
        if (count == 8 && selected.certified_closed_contour) {
            auto limited = selected_caps;
            limited.action_steps = 1;
            core::ResourcePlan denied_memory{UINT64_C(2000000000)};
            const auto denied = cf::solve_discrete_circulation_loads(
                input, write_loads, project_flow(), limited, selected, *geometry_work, shared,
                *local, *kernel, denied_memory);
            REQUIRE_FALSE(denied);
            CHECK(denied.error().code == core::FailureCode::ResourceExceeded);
        }
        const auto loads = take(cf::solve_discrete_circulation_loads(
            input, write_loads, std::move(flow), selected_caps, selected, *geometry_work, shared,
            *local, *kernel, memory));
        if (formal_product) {
#ifdef UVLAS_FINITE_WING_FORMAL_PRODUCT
            std::ifstream archive(UVLAS_FINITE_WING_COMMON_LIMIT);
            REQUIRE(archive);
            const auto frozen_limit = nlohmann::json::parse(archive);
            const auto& record = frozen_limit.at("records").at(UVLAS_FINITE_WING_INDEX);
            REQUIRE(record.at("dihedral") == wing.at("dihedral_degrees"));
            CHECK(loads.edge_count() == midpoint.size());
            check_formal_product(loads, grid, mesh, geometry, contact, transfers, sin[sample],
                                 record);
            std::cout << "Cumulative planned payload bytes " << memory.allocated_bytes() << '\n';
#endif
            return;
        }
        if (count == 8) {
            core::ResourcePlan replay_memory{UINT64_C(2000000000)};
            const auto retained_contact = take_global(contact.retain_checked(replay_memory));
            const auto contact_bytes = replay_memory.allocated_bytes();
            const cf::DiscreteNonpenetrationInput replay_input{
                retained_contact, rows, wake, input.history, trial, stamp, input.flow_side};
            const auto evidence = take(cf::verify_discrete_circulation_row_replay(
                replay_input, loads.primal(), selected_caps, selected, shared, *local, *kernel,
                replay_memory));
            CHECK(evidence.rows == count);
            CHECK(evidence.coefficients == count * count);
            CHECK(evidence.replay_row_bytes > 0);
            if (selected.certified_closed_contour)
                CHECK(evidence.reference_dense_bytes == 0);
            else
                CHECK(evidence.reference_dense_bytes == evidence.replay_row_bytes * count);
            CHECK(evidence.allocated_bytes == replay_memory.allocated_bytes() - contact_bytes);
            CHECK(evidence.used_steps > 0);
            std::cout << "Native coefficient replay, row bytes=" << evidence.replay_row_bytes
                      << ", dense bytes=" << evidence.reference_dense_bytes
                      << ", cumulative replay payload bytes=" << evidence.allocated_bytes
                      << ", used steps=" << evidence.used_steps << '\n';

            core::ResourcePlan insufficient{0};
            const auto denied = cf::verify_discrete_circulation_row_replay(
                input, loads.primal(), selected_caps, selected, shared, *local, *kernel,
                insufficient);
            if (denied)
                FAIL("The row replay must admit its payload before assembly");
            CHECK(denied.error().code == core::FailureCode::ResourceExceeded);
            CHECK(insufficient.allocated_bytes() == 0);
            auto changed = selected;
            changed.rhs_generation = core::RhsGeneration{9};
            const auto stale = cf::verify_discrete_circulation_row_replay(
                input, loads.primal(), selected_caps, changed, shared, *local, *kernel,
                insufficient);
            if (stale)
                FAIL("A different RHS generation must not qualify a replay");
            CHECK(stale.error().code == core::FailureCode::StaleGeneration);
            CHECK(insufficient.allocated_bytes() == 0);
            if (selected.certified_closed_contour) {
                auto parallel = selected;
                parallel.contour_workers = 4;
                parallel.solver_limits.dense_residual_workers = 4;
                const auto repeated = take(cf::solve_discrete_circulation_loads(
                    input, write_loads, project_flow(), selected_caps, parallel, *geometry_work,
                    shared, *local, *kernel, memory));
                const auto same_bits = [](double x, double y) {
                    CHECK(std::bit_cast<std::uint64_t>(x) == std::bit_cast<std::uint64_t>(y));
                };
                const auto equal_words = [&](const a::ScalarWords& cx, const a::ScalarWords& cy) {
                    REQUIRE(cx.index() == cy.index());
                    if (const auto* single = std::get_if<a::SingleWord>(&cx)) {
                        same_bits(single->value, std::get<a::SingleWord>(cy).value);
                    } else {
                        const auto x = std::get<a::DoubleWord>(cx);
                        const auto y = std::get<a::DoubleWord>(cy);
                        same_bits(x.high, y.high);
                        same_bits(x.low, y.low);
                    }
                };
                const auto equal_enclosure = [&](const a::ScaledEnclosure& x,
                                                 const a::ScaledEnclosure& y) {
                    CHECK(x.exponent() == y.exponent());
                    same_bits(x.radius(), y.radius());
                    equal_words(x.center(), y.center());
                };
                const auto equal_vector = [&](const cf::DiscreteLoadVector& x,
                                              const cf::DiscreteLoadVector& y) {
                    for (std::size_t axis = 0; axis < x.size(); ++axis)
                        equal_enclosure(x[axis], y[axis]);
                };
                CHECK(repeated.used_steps() == loads.used_steps());
                for (std::size_t i = 0; i < count; ++i) {
                    equal_words(words(loads.primal().solution(), i),
                                words(repeated.primal().solution(), i));
                    const auto xx = take(loads.primal().row(i));
                    const auto yy = take(repeated.primal().row(i));
                    CHECK(xx.face == yy.face);
                    equal_enclosure(xx.original_residual, yy.original_residual);
                }
                REQUIRE(repeated.edge_count() == loads.edge_count());
                for (std::size_t i = 0; i < loads.edge_count(); ++i) {
                    const auto x = take(loads.edge(i));
                    const auto y = take(repeated.edge(i));
                    CHECK(x.edge == y.edge);
                    CHECK(x.excluded_bound_occurrences == y.excluded_bound_occurrences);
                    equal_vector(x.force, y.force);
                }
                REQUIRE(repeated.point_count() == loads.point_count());
                for (std::size_t i = 0; i < loads.point_count(); ++i) {
                    const auto x = take(loads.point(i));
                    const auto y = take(repeated.point(i));
                    CHECK(x.point == y.point);
                    equal_vector(x.force, y.force);
                }
                equal_vector(loads.total_force(), repeated.total_force());
                equal_vector(loads.total_moment(), repeated.total_moment());
                auto too_tight = selected;
                too_tight.matrix_absolute = 0;
                auto unconsumed = project_flow();
                const auto unprojectable = cf::solve_discrete_circulation(
                    input, std::move(unconsumed), caps, too_tight, shared, *local, *kernel, memory);
                REQUIRE_FALSE(unprojectable);
                CHECK(unprojectable.error().code == core::FailureCode::PrecisionUnresolved);
                // NOLINTNEXTLINE(bugprone-use-after-move) A failed solve must not consume the flow.
                CHECK(unconsumed.source().jet.active());
            }
        }
        const auto tolerance = rational_word(2e-10);
        for (std::size_t k = 0; k < count; ++k) {
            const auto held = words(loads.primal().solution(), k);
            const Rat expected = sin[sample] * rational_word(gamma_prime->at(k).get<double>());
            CHECK(abs(exact_word(held) - expected) <= tolerance);
            check_within(unpack(take(loads.primal().row(k)).original_residual), 0, tolerance);
            rows[k].circulation = held;
        }
        // Recheck the original sheet with independent, explicitly supplied
        // downstream traces and the same retained two-word circulation.
        std::vector<cf::FrozenWakeTrace> traces;
        traces.reserve(ny);
        for (std::size_t k = 0; k < ny; ++k) {
            const auto& gamma = rows[count - ny + k].circulation;
            const auto value =
                transfers[k].coefficient == 1
                    ? gamma
                    : std::visit(
                          [](const auto& word) -> a::ScalarWords {
                              if constexpr (std::is_same_v<std::decay_t<decltype(word)>,
                                                           a::SingleWord>)
                                  return a::SingleWord{-word.value};
                              else
                                  return a::DoubleWord{-word.high, -word.low};
                          },
                          gamma);
            traces.push_back({a::P1StrengthId{700 + k}, value});
        }
        const cf::DiscreteNonpenetrationInput check_input{contact,
                                                          rows,
                                                          frozen_wake,
                                                          {core::AcceptedStateId{2}, traces},
                                                          trial,
                                                          stamp,
                                                          p::RequestedSide::Automatic};
        const auto closure = take(cf::evaluate_discrete_nonpenetration(
            check_input, project_flow(), caps, shared, *local, *kernel, memory));
        for (std::size_t k = 0; k < count; ++k)
            check_within(unpack(take(closure.row(k)).residual), 0, tolerance);

        REQUIRE(loads.edge_count() == midpoint.size());
        for (std::size_t k = 0; k < loads.edge_count(); ++k)
            CHECK(take(loads.edge(k)).excluded_bound_occurrences > 0);
        const auto force = unpack_vector(loads.total_force());
        near[sample] = add(scale(force[0], cos[sample]), scale(force[2], sin[sample]));
        far[sample] = trefftz_midpoint(loads.primal(), mesh, contact, grid.vertices,
                                       rational_word(1.2), count, ny);
    }
    REQUIRE(sin[0] == -sin[1]);
    REQUIRE(cos[0] == cos[1]);
    const Rat far_factor = 2 * sin[1] * sin[1];
    const Rat near_factor = far_factor * cos[1];
    REQUIRE(near_factor > 0);
    const auto near_coefficient = scale(add(near[0], near[1]), 1 / near_factor);
    const Rat far_coefficient = (far[0] + far[1]) / far_factor;
    const Rat finite_force_tolerance =
        rational_word(2e-11) * rational_word(1.2) * 20 * 20 * rational_word(8.4) / 2;
    check_within(near_coefficient, rational_word(oracle.at("near_drag_N_per_rad2").get<double>()),
                 4 * finite_force_tolerance / near_factor);
    check_within({far_coefficient, 0},
                 rational_word(oracle.at("trefftz_drag_N_per_rad2").get<double>()),
                 4 * finite_force_tolerance / far_factor);
    CHECK(far_coefficient > 0);
    std::cout << std::setprecision(16)
              << "Midpoint drag N/rad^2, near=" << near_coefficient.center.convert_to<double>()
              << ", near radius=" << near_coefficient.radius.convert_to<double>()
              << ", Trefftz=" << far_coefficient.convert_to<double>()
              << ", planned payload bytes=" << memory.allocated_bytes() << '\n';
    const Rat gap = abs(near_coefficient.center / far_coefficient - 1);
    const auto expected_gap = rational_word(oracle.at("relative_near_far_gap").get<double>());
    CHECK(abs(gap - expected_gap) <=
          (4 * finite_force_tolerance / near_factor + 4 * finite_force_tolerance / far_factor) /
              far_coefficient);
}
#endif

void verify_selected_wing() {
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
    static_assert(series >= 0 && series <= 3);
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
#ifdef UVLAS_FINITE_WING_REFINEMENT_FIXTURE
    } else if constexpr (series == 3) {
        std::ifstream polar_file{UVLAS_FINITE_WING_POLARS_FIXTURE};
        std::ifstream reference_file{UVLAS_FINITE_WING_REFINEMENT_FIXTURE};
        REQUIRE(polar_file);
        REQUIRE(reference_file);
        const auto polars = nlohmann::json::parse(polar_file);
        const auto reference = nlohmann::json::parse(reference_file);
        REQUIRE(reference.at("id") == "finite-wing-drag-refinement-v1");
#ifdef UVLAS_FINITE_WING_FORMAL_GEOMETRY
        REQUIRE(reference.at("scope") == "geometry_only");
#endif
        REQUIRE(reference.at("closure_m") == 1024.);
        REQUIRE(reference.at("core_m") == 1.114285714285714e-5);
        const auto& polar = polars.at("cases").at(wing_index);
        const auto& fixture = reference.at("cases").at(
            2 * static_cast<std::size_t>(UVLAS_FINITE_WING_REFINEMENT_LEVEL) + wing_index);
        REQUIRE(fixture.at("dihedral_degrees") == wing.at("dihedral_degrees"));
        run_refined_drag(wing, polar, fixture, reference.at("closure_m").get<double>(),
                         reference.at("core_m").get<double>());
#endif
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
} // namespace

#ifdef UVLAS_FORMAL_HARMONIC_FIXTURE
namespace {

std::vector<double> formal_harmonic_ages(std::span<const Point> collocation,
                                         std::span<const Point> trailing, double speed,
                                         double length, double omega) {
    double nearest = std::numeric_limits<double>::infinity();
    for (std::size_t i = 1; i < trailing.size(); ++i) {
        const auto& first = trailing[i - 1];
        const Point edge{trailing[i][0] - first[0], trailing[i][1] - first[1],
                         trailing[i][2] - first[2]};
        const double length_squared =
            std::inner_product(edge.begin(), edge.end(), edge.begin(), 0.);
        REQUIRE(length_squared > 0.);
        for (const auto& point : collocation) {
            const Point offset{point[0] - first[0], point[1] - first[1], point[2] - first[2]};
            const double t = std::clamp(
                std::inner_product(offset.begin(), offset.end(), edge.begin(), 0.) / length_squared,
                0., 1.);
            const double dx = offset[0] - t * edge[0];
            const double dy = offset[1] - t * edge[1];
            const double dz = offset[2] - t * edge[2];
            nearest = std::min(nearest, dx * dx + dy * dy + dz * dz);
        }
    }
    REQUIRE(nearest > 0.);
    auto step = std::min(.0005, std::sqrt(nearest) / (2 * speed));
    // The quadrature derives its final endpoint from the retained wake geometry.
    // Keep every supplied start strictly below the outward lower bound on that endpoint.
    const auto end_age_lower = std::nextafter(length / speed, 0.);
    const auto phase_limit = std::numbers::pi / omega;
    std::vector<double> starts{0.};
    for (auto next = std::min(step, phase_limit); next < end_age_lower;
         next = starts.back() + std::min(step, phase_limit)) {
        REQUIRE(next > starts.back());
        starts.push_back(next);
        step *= 1.5;
    }
    return starts;
}

void verify_formal_harmonic(double k, bool heave) {
    using Json = nlohmann::json;
    CAPTURE(k, heave);
    std::ifstream source{UVLAS_FORMAL_HARMONIC_FIXTURE};
    REQUIRE(source);
    const auto frozen = Json::parse(source);
    REQUIRE(frozen.at("id") == "finite-wing-unsteady-v1");
    REQUIRE(frozen.at("observations").size() == 12);
    REQUIRE(frozen.at("quadrature_order") == 16);
    const auto chord = k == .8 ? std::size_t{256} : std::size_t{128};
    constexpr std::size_t span = 192, width = span + 1;
    const auto count = chord * span;
    const auto order = static_cast<core::Index>(count);
    const auto& mesh_input = frozen.at("meshes").at(chord == 256 ? 1 : 0);
    REQUIRE(mesh_input.at("chord_panels") == chord);
    REQUIRE(mesh_input.at("span_panels") == span);
    auto vertices = mesh_input.at("vertices_m").get<std::vector<Point>>();
    const auto collocations = mesh_input.at("collocations_m").get<std::vector<Point>>();
    REQUIRE(vertices.size() == (chord + 1) * width);
    REQUIRE(collocations.size() == count);
    const auto wake_length = frozen.at("wake_length_m").get<double>();
    const auto trailing = std::vector<Point>(vertices.end() - width, vertices.end());
    for (auto tip : trailing) {
        tip[0] += wake_length;
        vertices.push_back(tip);
    }
    const FiniteWingMesh mesh{vertices, chord, span};
    std::vector<double> low(mesh.xyz.size());
    for (std::size_t i = 0; i < width; ++i) {
        const auto x = trailing[i][0];
        REQUIRE(std::abs(x) < wake_length);
        const auto tail = 3 * ((chord + 1) * width + i);
        low[tail] = x - (mesh.xyz[tail] - wake_length);
    }
    constexpr auto gib = UINT64_C(1) << 30U;
    core::ResourcePlan plan{24 * gib};
    const auto radius = k == .8 ? 1e-6 : 1e-5;
    const a::DoublePositions positions{mesh.xyz, low};
    std::vector<a::BandId> bands;
    bands.reserve(span);
    for (std::size_t i = 0; i < span; ++i)
        bands.emplace_back(50 + i);
    auto caps = test_geometry::limits();
    caps.counts = {mesh.ids.size(), mesh.edges.size(),     mesh.faces.size(),
                   bands.size(),    2 * mesh.faces.size(), mesh.quads.size()};
    caps.max_steps = 40000000;
    caps.predicate.max_steps = 1000000;
    auto geometry_work = test_geometry::workspace(plan);
    const Number instant{5, 4};
    const auto geometry =
        test_geometry::take(a::prepare_owned_geometry({mesh.ids,
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
                                                       instant.view()},
                                                      caps, *geometry_work, plan));
    const a::GlobalLimits global_caps{std::numeric_limits<std::uint64_t>::max() / 4,
                                      ls::EstimateAdmission{65536}};
    auto global = take(a::GlobalWorkspace::create(plan, global_caps.metadata));
    const auto topology = take_global(a::build_global_topology(
        geometry, {mesh.roles, mesh.owners, mesh.trailing}, global_caps, *global, plan));
    const auto contact =
        take_global(a::check_global_geometry(topology, geometry, global_caps, *global, plan));
    REQUIRE(topology.trailing_edges().size() == span);

    const auto speed = frozen.at("inflow_m_s").get<Point>();
    REQUIRE((speed == Point{20., 0., 0.}));
    const auto pivot = frozen.at("moment_origin_m").get<Point>();
    const auto omega = 2 * k * speed[0] / frozen.at("reduced_frequency_chord_m").get<double>();
    const auto ages = formal_harmonic_ages(collocations, trailing, speed[0], wake_length, omega);
    std::vector<a::ResolvedQuery> queries;
    std::vector<cf::DiscreteCollocation> rows;
    std::vector<cf::HarmonicCirculationRhs> rhs;
    std::vector<cf::HarmonicReflectionPair> pairs;
    queries.reserve(count);
    rows.reserve(count);
    rhs.reserve(count);
    pairs.reserve(count / 2);
    Rat rhs_squared{0};
    for (std::size_t i = 0; i < count; ++i) {
        const auto& x = collocations[i];
        queries.push_back(test_geometry::take(
            a::resolve_query(geometry, {a::QueryId{i + 1}, a::FreePoint{a::FreeSingle{x}}},
                             geometry.context(), caps, *geometry_work, plan)));
        rows.emplace_back(mesh.quads[i].facets[0], queries.back(),
                          cf::ExplicitCollocationTangents{xyz({1., 0., 0.}), xyz({0., 1., 0.})},
                          xyz({}), a::SingleWord{0}, a::SingleWord{0});
        const double real = heave ? 0. : -speed[0];
        const double imag = heave ? omega : -omega * (x[0] - pivot[0]);
        rhs.push_back({mesh.faces[i].id, a::SingleWord{real}, a::SingleWord{imag}});
        const auto r = rational_word(real), im = rational_word(imag);
        rhs_squared += r * r + im * im;
        if (i % span >= span / 2)
            pairs.push_back(
                {mesh.faces[i].id, mesh.faces[i / span * span + span - 1 - i % span].id});
    }
    std::vector<cf::DiscreteWakePanel> wake;
    wake.reserve(span);
    for (std::size_t i = 0; i < span; ++i)
        wake.emplace_back(mesh.quads[count + i].facets[0], cf::AttachedWakeFront{},
                          cf::SteadyWakeClosure{});
    constexpr a::P1ActionStamp trial{815, 6};
    const cf::DiscreteNonpenetrationInput input{contact,
                                                rows,
                                                wake,
                                                {core::AcceptedStateId{2}, {}},
                                                trial,
                                                stamp,
                                                p::RequestedSide::Automatic};
    const auto relative = frozen.at("relative_l2_residual_budget").get<double>();
    REQUIRE(relative == 1e-11);
    const std::vector<double> scales(count, 20.);
    const std::vector<ls::UnitDimensions> gamma_units(count, {2, 0, -1, 0, 0, 0, 0});
    const std::vector<ls::UnitDimensions> speed_units(count, {1, 0, -1, 0, 0, 0, 0});
    cf::DiscreteCirculationRequest request{
        {scales, scales, gamma_units, speed_units, context(geometry.context())},
        {order, 1, span, 5, ls::EstimateAdmission{14 * gib}},
        {ls::Norm::L2, 1e-12, 0., 1e-12, core::ObservationId{20}, core::ObservationId{21}},
        core::RhsGeneration{3},
        core::NumericalGeneration{4},
        6e-9,
        1e-13,
        .99 * relative * std::sqrt(rhs_squared.convert_to<double>()) / 20.,
        0.};
    request.certified_closed_contour = true;
    request.contour_workers = 4;
    request.solver_limits.dense_residual_workers = 4;
    const cf::DiscreteNonpenetrationLimits action_limits{
        {core::Precision::TwoWord, 1e-12},
        {core::Precision::TwoWord, 1e-12, 1e-12},
        {caps, {ls::EstimateAdmission{65536}, 1000000000}, UINT64_C(150000000000)},
        UINT64_C(150000000000),
        UINT64_C(70000000000000)};
    const cf::StraightWakeHarmonicInput phase{xyz(speed), ages, omega, 16 * ages.size(), pairs};
    const cf::HarmonicLoadInput loading{a::SingleWord{frozen.at("density_kg_m3").get<double>()},
                                        xyz(pivot)};
    const auto solved = take(cf::solve_discrete_harmonic_loads(input, phase, rhs, loading,
                                                               action_limits, request, plan));
    REQUIRE(solved.loads);
    REQUIRE(solved.residual.size() == count);
    REQUIRE(solved.physical_values().shape().rows == order);
    REQUIRE(solved.solution.values().shape().rows == order / 2);
    Rat squared{0};
    for (std::size_t i = 0; i < count; ++i) {
        const auto& row = solved.residual.values()[i];
        REQUIRE(row);
        CHECK(row->face == rhs[i].face);
        for (const auto& part : {unpack(row->real), unpack(row->imaginary)}) {
            const Rat upper = abs(part.center) + part.radius;
            squared += upper * upper;
        }
    }
    const auto tolerance = rational_word(relative);
    CHECK(squared < tolerance * tolerance * rhs_squared);
    std::cout << "Native relative L2 upper "
              << std::sqrt((squared / rhs_squared).convert_to<double>())
              << ", cumulative planned bytes " << plan.allocated_bytes() << '\n';

    std::size_t observations = 0;
    for (const auto& row : frozen.at("observations")) {
        if (row.at("k") != k || row.at("motion") != (heave ? "heave" : "pitch"))
            continue;
        ++observations;
        const bool force = row.at("observable") == "Fz";
        REQUIRE((force || row.at("observable") == "My"));
        REQUIRE(row.at("chord_panels") == chord);
        REQUIRE(row.at("span_panels") == span);
        REQUIRE(row.at("core_m") == radius);
        const auto real =
            unpack(force ? solved.loads->force_real[2] : solved.loads->moment_real[1]);
        const auto imag =
            unpack(force ? solved.loads->force_imaginary[2] : solved.loads->moment_imaginary[1]);
        const auto expected = row.at("reference_limit").get<std::array<double, 2>>();
        const auto reference_source = [&](double word) {
            const Rat center = rational_word(word);
            const auto infinity = std::numeric_limits<double>::infinity();
            const Rat upper = rational_word(std::nextafter(word, infinity)) - center;
            const Rat lower = center - rational_word(std::nextafter(word, -infinity));
            return Enclosure{center, std::max(lower, upper)};
        };
        const auto r = reference_source(expected[0]), im = reference_source(expected[1]);
        const Rat dr = abs(real.center - r.center) + real.radius + r.radius;
        const Rat di = abs(imag.center - im.center) + imag.radius + im.radius;
        const Rat reference_real =
            abs(r.center) > r.radius ? Rat{abs(r.center) - r.radius} : Rat{0};
        const Rat reference_imag =
            abs(im.center) > im.radius ? Rat{abs(im.center) - im.radius} : Rat{0};
        const Rat reference_norm_squared =
            reference_real * reference_real + reference_imag * reference_imag;
        REQUIRE(reference_norm_squared > 0);
        const Rat uncertainty =
            rational_word(std::nextafter(row.at("reference_relative_uncertainty").get<double>(),
                                         std::numeric_limits<double>::infinity()));
        CHECK(uncertainty < Rat{1} / 400);
        REQUIRE(frozen.at("combined_observation_budget") == .01);
        const Rat remaining = Rat{1} / 100 - uncertainty;
        REQUIRE(remaining > 0);
        CHECK(dr * dr + di * di < remaining * remaining * reference_norm_squared);
        std::cout << std::setprecision(16) << row.at("observable") << " = "
                  << real.center.convert_to<double>() << ", " << imag.center.convert_to<double>()
                  << "; combined relative upper "
                  << std::sqrt(
                         ((dr * dr + di * di) / reference_norm_squared).convert_to<double>()) +
                         uncertainty.convert_to<double>()
                  << '\n';
    }
    CHECK(observations == 2);
}

} // namespace

TEST_CASE("Formal harmonic wing k=0.05", "[.formal-harmonic][k005][V-AERO][verification]") {
    verify_formal_harmonic(.05, true);
    verify_formal_harmonic(.05, false);
}

TEST_CASE("Formal harmonic wing k=0.2", "[.formal-harmonic][k020][V-AERO][verification]") {
    verify_formal_harmonic(.2, true);
    verify_formal_harmonic(.2, false);
}

TEST_CASE("Formal harmonic wing k=0.8", "[.formal-harmonic][k080][V-AERO][verification]") {
    verify_formal_harmonic(.8, true);
    verify_formal_harmonic(.8, false);
}
#endif

#ifdef UVLAS_FINITE_WING_FORMAL_PRODUCT
TEST_CASE("A formal finite-wing grid qualifies steady circulation and induced drag",
          "[V-AERO][AE-02][UW-D03][UW-D05][N4-06][verification]") {
    verify_selected_wing();
}
#elif defined(UVLAS_FINITE_WING_FORMAL_GEOMETRY)
TEST_CASE("A formal finite-wing grid qualifies geometry and single-query forecasts",
          "[V-AERO][AE-02][UW-D03][UW-D05][N4-06][verification]") {
    verify_selected_wing();
}
#else
TEST_CASE("A frozen finite wing reproduces its complete steady polar and loads",
          "[V-AERO][AE-02][UW-D03][UW-D05][N4-06][verification]") {
    verify_selected_wing();
}
#endif
