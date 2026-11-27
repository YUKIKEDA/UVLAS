#include "case_pg_shell_material.hpp"

#include "io/si_case_detail.hpp"
#include "structure/pg_reference_arithmetic.hpp"

namespace uvlas::coupling::detail {
namespace {
namespace pg = structure;
using core::FailureCode;
using pg::detail::partition_entry;
using pg::detail::partition_fail;
using pg::detail::partition_index;
using pg::detail::partition_take;
using pg::detail::ReferenceArithmetic;
constexpr auto operation = core::Operation::PrepareModel;
constexpr io::SchemaDimension pressure{{-1, 1, -2, 0}};
constexpr io::SchemaDimension density{{-3, 1, 0, 0}};
constexpr io::SchemaDimension dimensionless{{0, 0, 0, 0}};

template <class Field>
pg::ExactRational quantity(const io::SiCase& source, Field field, io::SchemaDimension dimension,
                           ReferenceArithmetic& math) {
    math.charge();
    const auto recipe = partition_take(source.quantity(field)).recipe;
    if (recipe.dimensions.exponents != dimension.exponents || recipe.pi_power != 0)
        partition_fail(FailureCode::InvalidInput);
    return partition_take(io::materialize_rational(recipe, math.leaf_limits, math.plan));
}

template <class Field>
std::array<pg::ExactRational, 3> triple(const io::SiCase& source, Field field,
                                        io::SchemaDimension dimension, ReferenceArithmetic& math) {
    std::array<pg::ExactRational, 3> values;
    for (std::size_t i = 0; i < values.size(); ++i)
        values[i] = quantity(source, partition_take(field.at(partition_index(i))), dimension, math);
    return values;
}

// Worst-case owned rational slots per reusable law (including both Mandel coefficients).
constexpr std::uint64_t law_rationals = 2 * 36 + 2;

pg::PgShellElasticLaw prepare_law(const io::SiCase& source, const auto& record,
                                  std::uint64_t other_rationals, ReferenceArithmetic& math) {
    using Kind = std::decay_t<decltype(record)>;
    const auto rho = quantity(source, record.density, density, math);
    const auto qualify = [&](const pg::PgElasticity& input) {
        auto remaining = math.limits;
        remaining.maximum_scalar_operations -= math.scalar_operations;
        math.reserve_nested_rationals(remaining.algebra, other_rationals);
        auto law = partition_take(pg::PgShellElasticLaw::prepare(input, rho, remaining, math.plan));
        math.charge(law.call_info().scalar_operations);
        return law;
    };
    if constexpr (std::is_same_v<Kind, io::schema::MaterialIsotropic>) {
        const auto young = quantity(source, record.young, pressure, math);
        const auto poisson = quantity(source, record.poisson, dimensionless, math);
        return qualify(pg::PgIsotropicElasticity{young, poisson});
    } else if constexpr (std::is_same_v<Kind, io::schema::MaterialOrthotropic>) {
        const auto young = triple(source, record.young, pressure, math);
        const auto poisson = triple(source, record.poisson, dimensionless, math);
        const auto shear = triple(source, record.shear, pressure, math);
        return qualify(pg::PgOrthotropicElasticity{young, poisson, shear});
    } else {
        std::array<pg::ExactRational, 36> stiffness;
        for (core::Index i = 0; i < 6; ++i) {
            const auto row = partition_take(record.stiffness.at(i));
            for (core::Index j = 0; j < 6; ++j) {
                stiffness[static_cast<std::size_t>(6 * i + j)] =
                    quantity(source, partition_take(row.at(j)), pressure, math);
            }
        }
        const auto tag = record.components.selected();
        const auto convention = tag == decltype(tag)::ItemEngineering
                                    ? pg::PgMaterialComponents::Engineering
                                    : pg::PgMaterialComponents::Mandel;
        return qualify(pg::PgEffectiveElasticity{stiffness, convention});
    }
}
} // namespace

CasePgShellMaterials::CasePgShellMaterials(io::SchemaIdentity source, core::Index shell,
                                           core::PlannedBuffer<CasePgShellLayer>&& layers,
                                           Laws&& laws, std::size_t law_count) noexcept
    : source_(source),
      shell_(shell),
      layers_(std::move(layers)),
      laws_(std::move(laws)),
      law_count_(law_count) {}

core::Result<CasePgShellMaterials>
CasePgShellMaterials::prepare(const io::SiCase& source, std::uint64_t shell_id,
                              const pg::PgReferenceLimits& limits,
                              core::ResourcePlan& memory) noexcept {
    return partition_entry([&]() -> core::Result<CasePgShellMaterials> {
        ReferenceArithmetic math{limits, memory};
        math.cancel();
        const auto& typed = partition_take(source.source_case()).get();
        const auto identity = partition_take(typed.identity());
        const auto model = partition_take(typed.root()).model.record();
        core::Index shell_ordinal = 0;
        for (; shell_ordinal < partition_index(model.shells.size()); ++shell_ordinal) {
            math.visit(1);
            if (partition_take(model.shells.at(shell_ordinal)).record().id.unsigned_value() ==
                shell_id)
                break;
        }
        const auto shell = partition_take(model.shells.at(shell_ordinal)).record();
        const auto count = shell.layers.size();
        if (count == 0)
            partition_fail(FailureCode::InvalidInput);
        const auto bounds_rationals = pg::detail::partition_mul(count, 2);
        const auto all_laws =
            pg::detail::partition_mul(law_rationals, std::min(count, model.materials.size()));
        math.admit(
            pg::detail::partition_add(pg::detail::reference_live_rationals,
                                      pg::detail::partition_add(bounds_rationals, all_laws)));
        auto layers =
            partition_take(memory.allocate<CasePgShellLayer>(partition_index(count), operation));
        auto laws = partition_take(memory.allocate<std::optional<pg::PgShellElasticLaw>>(
            partition_index(std::min(count, model.materials.size())), operation));
        struct Index {
            std::uint64_t id{};
            core::Index ordinal{};
        };
        auto by_material = partition_take(
            memory.allocate<Index>(partition_index(model.materials.size()), operation));
        auto by_layer = partition_take(memory.allocate<Index>(partition_index(count), operation));
        const auto less = [&](std::uint64_t a, std::uint64_t b) {
            math.visit(1);
            return a < b;
        };
        for (std::size_t i = 0; i < by_material.size(); ++i) {
            math.visit(1);
            const auto ordinal = partition_index(i);
            const auto material = partition_take(model.materials.at(ordinal)).record();
            by_material.values()[i] = {
                std::visit([](const auto& value) { return value.id.unsigned_value(); }, material),
                ordinal};
        }
        auto lower = math.constant("-1");
        for (std::size_t i = 0; i < layers.size(); ++i) {
            math.visit(1);
            const auto ordinal = partition_index(i);
            const auto layer = partition_take(shell.layers.at(ordinal)).record();
            auto lo = quantity(source, partition_take(layer.interval.at(0)), dimensionless, math);
            auto hi = quantity(source, partition_take(layer.interval.at(1)), dimensionless, math);
            if (!(lo == lower) || math.compare(hi, lo) <= 0)
                partition_fail(FailureCode::InvalidInput);
            lower = math.copy(hi);
            const auto material_id = layer.material.unsigned_value();
            layers.values()[i] = {core::LayerId{layer.id.unsigned_value()}, material_id,
                                  std::move(lo), std::move(hi), 0};
            by_layer.values()[i] = {material_id, ordinal};
        }
        if (!(lower == math.constant("1")))
            partition_fail(FailureCode::InvalidInput);

        const auto index = by_material.values();
        std::ranges::sort(index, less, &Index::id);
        std::ranges::sort(by_layer.values(), less, &Index::id);
        std::optional<std::uint64_t> previous;
        std::size_t prepared = 0;
        for (const auto& entry : by_layer.values()) {
            math.visit(1);
            if (entry.id != previous) {
                const auto found = std::ranges::lower_bound(index, entry.id, less, &Index::id);
                if (found == index.end() || found->id != entry.id || prepared >= laws.size())
                    partition_fail(FailureCode::InvalidInput);
                const auto record = partition_take(model.materials.at(found->ordinal)).record();
                const auto held = pg::detail::partition_add(
                    bounds_rationals, pg::detail::partition_mul(prepared, law_rationals));
                laws.values()[prepared++].emplace(std::visit(
                    [&](const auto& value) { return prepare_law(source, value, held, math); },
                    record));
                previous = entry.id;
            }
            layers.values()[static_cast<std::size_t>(entry.ordinal)].law = prepared - 1;
        }
        math.cancel();
        return CasePgShellMaterials{identity, shell_ordinal, std::move(layers), std::move(laws),
                                    prepared};
    });
}

core::Result<std::reference_wrapper<const pg::PgShellElasticLaw>>
CasePgShellMaterials::law_for(core::Index layer) const& noexcept {
    if (layer < 0 || static_cast<std::uint64_t>(layer) >= layers_.size())
        return std::unexpected(core::Failure{FailureCode::InvalidInput, operation});
    const auto index = layers_.values()[static_cast<std::size_t>(layer)].law;
    return std::cref(laws_.values()[index].value());
}

core::Result<io::schema::Layer>
CasePgShellMaterials::source_layer(const io::SiCase& source, core::Index layer) const noexcept {
    return partition_entry([&]() -> core::Result<io::schema::Layer> {
        static_cast<void>(partition_take(law_for(layer)));
        const auto& typed = partition_take(source.source_case()).get();
        if (!io::detail::same_identity(partition_take(typed.identity()), source_))
            partition_fail(FailureCode::StaleGeneration);
        const auto model = partition_take(typed.root()).model.record();
        const auto shell = partition_take(model.shells.at(shell_)).record();
        return partition_take(shell.layers.at(layer)).record();
    });
}
} // namespace uvlas::coupling::detail
