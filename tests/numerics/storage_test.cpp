#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <uvlas/numerics/storage.hpp>

namespace {
using namespace uvlas::core;
using namespace uvlas::numerics;
constexpr Operation op = Operation::PrepareModel;
static_assert(StorageScalar<double>);
static_assert(StorageScalar<std::complex<double>>);
static_assert(!StorageScalar<float>);
static_assert(!StorageScalar<long double>);
constexpr double quiet_nan = std::numeric_limits<double>::quiet_NaN();
template <class T>
void failure(const Result<T>& result, FailureCode code) {
    REQUIRE_FALSE(result);
    CHECK(result.error().code == code);
    CHECK(result.error().operation == op);
}
constexpr SiDimension length_unit{{1, 0, 0, 0, 0, 0, 0}};
constexpr SiDimension no_unit{{0, 0, 0, 0, 0, 0, 0}};
std::array<DofDescriptor, 3> descriptors() {
    return {{{DofId{91}, ComponentId{2}, CoordinateKind::Position, length_unit, 2.5},
             {DofId{11}, ComponentId{7}, CoordinateKind::AncfGradient, no_unit, 0.125},
             {DofId{55}, ComponentId{3}, CoordinateKind::BeamRotation, no_unit, 3.0}}};
}
std::array<DofAssociation, 4> owners() {
    return {{{2, MemberId{8}}, {0, PatchId{3}}, {1, PatchId{3}}, {0, PatchId{9}}}};
}
} // namespace

TEST_CASE("Dense copies preserve padding complex parts and ownership", "[unit][N4-02]") {
    ResourcePlan plan{4096};
    std::array<double, 6> source{1, 2, 99, 3, 4, 88};
    auto real = DenseReal::copy({2, 2, 3}, source, plan, op);
    REQUIRE(real);
    CHECK(real->shape().leading_dimension == 3);
    CHECK(real->values()[4] == 4);
    source[0] = 17;
    CHECK(real->values()[0] == 1);
    const std::array<std::complex<double>, 2> data{{{1, 2}, {3, -4}}};
    auto complex = DenseComplex::copy({2, 1, 2}, data, plan, op);
    REQUIRE(complex);
    CHECK(complex->values()[1] == std::complex<double>{3, -4});
    DenseReal moved{std::move(*real)};
    CHECK(real->shape().rows == 0);
    CHECK(real->values().empty());
    CHECK(moved.values()[0] == 1);
    moved.values()[0] = quiet_nan;
    failure(moved.validate(op), FailureCode::Nonfinite);
    STATIC_REQUIRE_FALSE(std::is_copy_constructible_v<DenseReal>);
}

TEST_CASE("Dense rejects invalid storage and accepts empty shapes", "[negative][N4-02]") {
    ResourcePlan plan{4096};
    const std::array<double, 2> data{1, quiet_nan};
    failure(DenseReal::copy({1, 1, 2}, data, plan, op), FailureCode::Nonfinite);
    failure(DenseReal::copy({2, 1, 1}, data, plan, op), FailureCode::InvalidInput);
    failure(DenseReal::copy({-1, 1, 1}, {}, plan, op), FailureCode::InvalidInput);
    failure(DenseReal::copy({1, 2, std::numeric_limits<Index>::max()}, {}, plan, op),
            FailureCode::ResourceExceeded);
    for (const auto value :
         {std::complex<double>{1, quiet_nan}, std::complex<double>{quiet_nan, 1}}) {
        const std::array<std::complex<double>, 1> complex{value};
        failure(DenseComplex::copy({1, 1, 1}, complex, plan, op), FailureCode::Nonfinite);
    }
    ResourcePlan empty_plan{0};
    for (const auto shape : {DenseShape{0, 0, 0}, DenseShape{0, 3, 7}, DenseShape{4, 0, 4}}) {
        const auto empty = DenseReal::copy(shape, {}, empty_plan, op);
        REQUIRE(empty);
        CHECK(empty->values().empty());
    }
    const std::array<double, 1> one{1};
    failure(DenseReal::copy({1, 1, 1}, one, empty_plan, op), FailureCode::ResourceExceeded);
}

TEST_CASE("Dense zero initialization charges only the final matrix", "[unit][negative][N4-02]") {
    constexpr DenseShape padded{2, 2, 3};
    constexpr auto bytes = 6 * sizeof(double);
    ResourcePlan short_plan{bytes - 1};
    failure(DenseReal::zeroed(padded, short_plan, op), FailureCode::ResourceExceeded);
    CHECK(short_plan.allocated_bytes() == 0);

    ResourcePlan exact{bytes};
    const auto initialized = DenseReal::zeroed(padded, exact, op);
    REQUIRE(initialized);
    CHECK(initialized->shape().leading_dimension == 3);
    REQUIRE(initialized->values().size() == 6);
    for (const double value : initialized->values()) {
        CHECK(value == 0);
        CHECK_FALSE(std::signbit(value));
    }
    CHECK(exact.allocated_bytes() == bytes);
    failure(DenseReal::zeroed({2, 1, 1}, exact, op), FailureCode::InvalidInput);
    const auto empty = DenseReal::zeroed({0, 3, 7}, exact, op);
    REQUIRE(empty);
    CHECK(empty->values().empty());
}

TEST_CASE("Sparse compression preserves compensated duplicate sums", "[unit][N4-02]") {
    ResourcePlan plan{8192};
    const std::array<Triplet, 6> input{
        {{2, 1, 8}, {0, 0, 1e16}, {1, 0, 3}, {0, 0, 1}, {0, 0, -1e16}, {2, 1, -8}}};
    const auto matrix = SparseReal::from_triplets(3, 2, input, plan, op);
    REQUIRE(matrix);
    CHECK(matrix->column_offsets().size() == 3);
    CHECK(matrix->column_offsets()[0] == 0);
    CHECK(matrix->column_offsets()[1] == 2);
    CHECK(matrix->column_offsets()[2] == 3);
    CHECK(matrix->row_indices()[0] == 0);
    CHECK(matrix->row_indices()[1] == 1);
    CHECK(matrix->row_indices()[2] == 2);
    // Exact integer sums, independent of the aggregation implementation.
    CHECK(matrix->values()[0] == 1);
    CHECK(matrix->values()[1] == 3);
    CHECK(matrix->values()[2] == 0);
}

TEST_CASE("Sparse CSC rejects malformed patterns and nonfinite aggregation", "[negative][N4-02]") {
    ResourcePlan plan{8192};
    const std::array<Index, 3> offsets{0, 2, 2};
    const std::array<Index, 2> indices{0, 1};
    const std::array<double, 2> values{2, 4};
    auto valid = SparseReal::copy_csc(2, 2, offsets, indices, values, plan, op);
    REQUIRE(valid);
    for (const auto bad : {std::array<Index, 3>{0, -1, 2}, std::array<Index, 3>{0, 3, 2},
                           std::array<Index, 3>{1, 2, 2}, std::array<Index, 3>{0, 1, 1}}) {
        failure(SparseReal::copy_csc(2, 2, bad, indices, values, plan, op),
                FailureCode::InvalidInput);
    }
    for (const auto bad : {std::array<Index, 2>{1, 0}, std::array<Index, 2>{0, 0},
                           std::array<Index, 2>{-1, 0}, std::array<Index, 2>{0, 2}}) {
        failure(SparseReal::copy_csc(2, 2, offsets, bad, values, plan, op),
                FailureCode::InvalidInput);
    }
    const std::array<double, 2> nonfinite{2, quiet_nan};
    failure(SparseReal::copy_csc(2, 2, offsets, indices, nonfinite, plan, op),
            FailureCode::Nonfinite);
    const std::array<Triplet, 2> overflow{
        {{0, 0, std::numeric_limits<double>::max()}, {0, 0, std::numeric_limits<double>::max()}}};
    failure(SparseReal::from_triplets(1, 1, overflow, plan, op), FailureCode::Nonfinite);
    const std::array<Triplet, 1> negative{{{-1, 0, 1}}};
    failure(SparseReal::from_triplets(1, 1, negative, plan, op), FailureCode::InvalidInput);
    failure(SparseReal::from_triplets(0, std::numeric_limits<Index>::max(), {}, plan, op),
            FailureCode::ResourceExceeded);
    ResourcePlan zero{0};
    const auto empty = SparseReal::from_triplets(0, 0, {}, zero, op);
    REQUIRE(empty);
    CHECK(empty->column_offsets().size() == 1);
    CHECK(empty->column_offsets()[0] == 0);
    CHECK(zero.allocated_bytes() == 0);
}

TEST_CASE("PointSoA owns matching coordinates and unique IDs", "[negative][N4-02]") {
    ResourcePlan plan{4096};
    std::array<double, 2> x{1, 2};
    const std::array<double, 2> y{3, 4}, z{5, 6};
    const std::array<PointId, 2> ids{PointId{0},
                                     PointId{std::numeric_limits<std::uint64_t>::max()}};
    auto points = PointSoA::copy(x, y, z, ids, plan, op);
    REQUIRE(points);
    x[0] = 9;
    CHECK(points->x()[0] == 1);
    CHECK(points->z()[1] == 6);
    REQUIRE(points->id(1, op));
    CHECK(points->id(1, op)->value() == std::numeric_limits<std::uint64_t>::max());
    failure(points->id(-1, op), FailureCode::InvalidInput);
    failure(points->id(2, op), FailureCode::InvalidInput);
    failure(PointSoA::copy(x, y, std::span(z).first(1), ids, plan, op), FailureCode::InvalidInput);
    const std::array<PointId, 2> duplicates{PointId{2}, PointId{2}};
    failure(PointSoA::copy(x, y, z, duplicates, plan, op), FailureCode::InvalidInput);
    for (int axis = 0; axis < 3; ++axis) {
        const std::array<double, 2> bad{quiet_nan, 2};
        failure(PointSoA::copy(axis == 0 ? std::span<const double>(bad) : x,
                               axis == 1 ? std::span<const double>(bad) : y,
                               axis == 2 ? std::span<const double>(bad) : z, ids, plan, op),
                FailureCode::Nonfinite);
    }
    ResourcePlan zero{0};
    REQUIRE(PointSoA::copy({}, {}, {}, {}, zero, op));
}

TEST_CASE("DofLayout owns metadata topology and coordinate selection", "[unit][N4-02]") {
    ResourcePlan plan{8192};
    auto dofs = descriptors();
    auto associations = owners();
    const std::array<Index, 2> independent{2, 0};
    auto layout = DofLayout::copy(dofs, independent, associations, plan, op);
    REQUIRE(layout);
    dofs[0].scale = 99;
    associations[0].full_index = 0;
    REQUIRE(layout->descriptor(0, op));
    CHECK(layout->descriptor(0, op)->scale == 2.5);
    CHECK(layout->descriptor(0, op)->id == DofId{91});
    CHECK(layout->descriptor(0, op)->unit.powers == length_unit.powers);
    CHECK(layout->descriptor(1, op)->kind == CoordinateKind::AncfGradient);
    CHECK(layout->descriptor(1, op)->id == DofId{11});
    CHECK(layout->descriptor(1, op)->component == ComponentId{7});
    CHECK(layout->descriptor(1, op)->unit.powers == no_unit.powers);
    CHECK(layout->descriptor(1, op)->scale == 0.125);
    CHECK(layout->descriptor(2, op)->kind == CoordinateKind::BeamRotation);
    CHECK(layout->descriptor(2, op)->id == DofId{55});
    CHECK(layout->descriptor(2, op)->component == ComponentId{3});
    CHECK(layout->descriptor(2, op)->scale == 3.0);
    CHECK(layout->independent_to_full(0, op) == Index{2});
    CHECK(layout->independent_to_full(1, op) == Index{0});
    REQUIRE(layout->full_to_independent(0, op));
    CHECK(*layout->full_to_independent(0, op) == std::optional<Index>{1});
    REQUIRE(layout->full_to_independent(1, op));
    CHECK_FALSE(layout->full_to_independent(1, op)->has_value());
    CHECK(layout->association_size() == 4);
    CHECK(layout->association(3, op)->full_index == 2);
    CHECK(std::holds_alternative<MemberId>(layout->association(3, op)->owner));
    CHECK(layout->association(0, op)->owner == DofOwner{PatchId{3}});
    CHECK(layout->association(1, op)->owner == DofOwner{PatchId{9}});
    CHECK(layout->association(2, op)->owner == DofOwner{PatchId{3}});
    CHECK(layout->association(3, op)->owner == DofOwner{MemberId{8}});
    failure(layout->descriptor(-1, op), FailureCode::InvalidInput);
    failure(layout->independent_to_full(2, op), FailureCode::InvalidInput);
    failure(layout->full_to_independent(3, op), FailureCode::InvalidInput);
    auto all_dependent = DofLayout::copy(descriptors(), {}, owners(), plan, op);
    REQUIRE(all_dependent);
    CHECK(all_dependent->independent_size() == 0);
    REQUIRE(all_dependent->full_to_independent(0, op));
    CHECK_FALSE(all_dependent->full_to_independent(0, op)->has_value());
}

TEST_CASE("DofLayout rejects invalid metadata and mappings", "[negative][N4-02]") {
    const auto full = descriptors();
    const auto associations = owners();
    for (const auto map :
         {std::array<Index, 2>{0, 0}, std::array<Index, 2>{-1, 1}, std::array<Index, 2>{1, 3}}) {
        ResourcePlan plan{8192};
        failure(DofLayout::copy(full, map, associations, plan, op), FailureCode::InvalidInput);
    }
    for (const double scale : {0.0, -1.0, quiet_nan, std::numeric_limits<double>::infinity()}) {
        auto bad = full;
        bad[1].scale = scale;
        ResourcePlan plan{8192};
        failure(DofLayout::copy(bad, {}, associations, plan, op),
                std::isfinite(scale) ? FailureCode::InvalidInput : FailureCode::Nonfinite);
    }
    auto duplicates = full;
    duplicates[1].id = duplicates[0].id;
    ResourcePlan plan{16384};
    failure(DofLayout::copy(duplicates, {}, associations, plan, op), FailureCode::InvalidInput);
    failure(DofLayout::copy(full, {}, std::span(associations).first(1), plan, op),
            FailureCode::InvalidInput);
    auto duplicated_owner = associations;
    duplicated_owner[3] = duplicated_owner[1];
    failure(DofLayout::copy(full, {}, duplicated_owner, plan, op), FailureCode::InvalidInput);
    auto invalid_owner = associations;
    invalid_owner[0].full_index = 3;
    failure(DofLayout::copy(full, {}, invalid_owner, plan, op), FailureCode::InvalidInput);
    auto invalid_kind = full;
    invalid_kind[0].kind = static_cast<CoordinateKind>(255);
    failure(DofLayout::copy(invalid_kind, {}, associations, plan, op), FailureCode::InvalidInput);
    ResourcePlan zero{0};
    REQUIRE(DofLayout::copy({}, {}, {}, zero, op));
}
