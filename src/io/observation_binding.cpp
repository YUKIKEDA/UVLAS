#include "observation_binding.hpp"

#include "source_digest.hpp"

#include <algorithm>
#include <type_traits>
#include <uvlas/io/si_case.hpp>

namespace uvlas::io::detail {
namespace {

using Code = core::FailureCode;
using Exact = structure::ExactRational;
constexpr auto operation = core::Operation::PrepareModel;
constexpr SchemaDimension length{{1, 0, 0, 0}}, time{{0, 0, 1, 0}};
constexpr SchemaDimension force{{1, 1, -2, 0}}, moment{{2, 1, -2, 0}}, unitless{};
constexpr std::array<std::string_view, 3> xyz{"x", "y", "z"};
constexpr std::array<std::string_view, 6> wrench{"Fx", "Fy", "Fz", "Mx", "My", "Mz"};

void require(bool condition, Code failure = Code::InvalidInput) {
    if (!condition)
        throw core::Failure{failure, operation};
}

template <class T>
T take(core::Result<T> value) {
    if (!value)
        throw value.error();
    return std::move(*value);
}

struct Binding {
    Binding(const ObservationFile& input, const Exact& time_value, const std::stop_token& cancel,
            SiUnits input_units, structure::ExactLimits arithmetic, std::uint64_t work_limit,
            core::ResourcePlan& resource, ObservationBindingDiagnostic where) noexcept
        : file(input),
          evaluation_time(time_value),
          stop(cancel),
          units(input_units),
          exact(arithmetic),
          work(work_limit),
          plan(resource),
          diagnostic(where) {}

    const ObservationFile& file;
    const Exact& evaluation_time;
    const std::stop_token& stop;
    SiUnits units;
    structure::ExactLimits exact;
    std::uint64_t work;
    core::ResourcePlan& plan;
    ObservationBindingDiagnostic diagnostic;
    std::size_t bound = 0;

    void step(std::uint64_t amount = 1) {
        require(!stop.stop_requested(), Code::Cancelled);
        require(amount <= work, Code::ResourceExceeded);
        work -= amount;
    }

    Exact converted(SchemaNumberView value, SchemaDimension dimensions) {
        step(exact.maximum_work_bits);
        const SiRecipeView recipe{
            SiDecimalView{value.sign, value.digits, value.exponent10}, {}, units, dimensions, 0};
        auto result = take(materialize_rational(recipe, exact, plan));
        if (auto* memory = std::get_if<structure::ExactEstimateAdmission>(&exact.external_memory)) {
            const auto used = result.preparation().estimated_external_payload_bytes;
            require(used <= memory->bytes, Code::ResourceExceeded);
            memory->bytes -= used;
        }
        return result;
    }

    const ObservationFileRecord& find(core::ObservationId id) {
        const auto records = file.records.values();
        const auto item = std::ranges::lower_bound(records, id.value(), {}, [&](const auto& x) {
            step();
            return x.id.value();
        });
        require(item != records.end() && item->id == id);
        return *item;
    }

    template <class ComponentList>
    void match_vector(const ComponentList& fields, const std::array<Exact, 3>& expected,
                      SchemaDimension dimensions) {
        for (core::Index i = 0; i < 3; ++i)
            require(converted(take(fields.at(i)).number(), dimensions) ==
                    expected[static_cast<std::size_t>(i)]);
    }

    template <class Observation>
    void match_descriptor(const Observation& observation, const ObservationFileRecord& row) {
        const auto basis = observation.basis.record();
        const auto* world = std::get_if<schema::BasisWorldMaterial>(&basis);
        require(world && world->kind.text() == "world", Code::Unsupported);
        if constexpr (std::is_same_v<Observation, schema::ObservationPoint>) {
            require(observation.quantity.text() == "displacement", Code::Unsupported);
            const auto point = observation.point.record();
            const auto target = point.target.record();
            require(target.kind.text() == "shell", Code::Unsupported);
            const auto* expected = std::get_if<ObservationFilePoint>(&row.descriptor);
            require(expected && expected->shell == target.id.unsigned_value());
            const auto side = static_cast<ObservationPointSide>(point.side.selected());
            require(side == expected->side);
            match_vector(point.coordinates, expected->coordinate, unitless);
        } else {
            const auto* expected = std::get_if<ObservationFileReaction>(&row.descriptor);
            require(expected && expected->support == observation.support.unsigned_value());
            match_vector(observation.origin, expected->origin_metres, length);
        }
    }

    void match_values(const schema::ComparisonReferenceFile& declaration,
                      const ObservationFileRecord& row, bool point) {
        const auto names = point ? std::span<const std::string_view>{xyz}
                                 : std::span<const std::string_view>{wrench};
        const auto count = declaration.columns.size();
        require(count <= names.size() && count == declaration.uncertainty.size());
        std::array<bool, 6> selected{};
        for (std::uint64_t i = 0; i < count; ++i) {
            step();
            const auto index = static_cast<core::Index>(i);
            const auto name = take(declaration.columns.at(index)).text();
            step(name.size());
            const auto item = std::ranges::find(names, name);
            require(item != names.end());
            const auto component = static_cast<std::size_t>(item - names.begin());
            const auto& measured = row.components[component];
            if (!measured || selected[component])
                throw core::Failure{Code::InvalidInput, operation};
            const auto uncertainty = take(declaration.uncertainty.at(index)).number();
            const auto dimensions = point ? length : component < 3 ? force : moment;
            require(converted(uncertainty, dimensions) == measured->uncertainty);
            selected[component] = true;
        }
        for (std::size_t i = 0; i < row.components.size(); ++i)
            require(selected[i] == row.components[i].has_value());
    }

    template <class Observation>
    void match(const Observation& observation, core::Index index) {
        step();
        const auto comparison = observation.comparison.record();
        const auto reference = comparison.reference.record();
        const auto* declaration = std::get_if<schema::ComparisonReferenceFile>(&reference);
        if (!declaration)
            return;
        const auto source = declaration->source.record();
        if (source.source_id.unsigned_value() != file.source.source.value())
            return;

        diagnostic.observation = core::ObservationId{observation.id.unsigned_value()};
        diagnostic.case_index = index;
        require(source.role.text() == "observation-reference");
        require(decode_source_digest(source.sha256.text()) == file.source_sha256);
        if constexpr (std::is_same_v<Observation, schema::ObservationPoint> ||
                      std::is_same_v<Observation, schema::ObservationReaction>) {
            require(comparison.statistic.text() == "instantaneous" &&
                        comparison.side.text() == "right",
                    Code::Unsupported);
            const auto& row = find(*diagnostic.observation);
            diagnostic.input_node = row.input_node;
            match_descriptor(observation, row);
            const auto left = converted(take(comparison.window.at(0)).number(), time);
            const auto right = converted(take(comparison.window.at(1)).number(), time);
            require(left == right, Code::Unsupported);
            require(left == evaluation_time && row.time_seconds == evaluation_time);
            match_values(*declaration, row, std::is_same_v<Observation, schema::ObservationPoint>);
            ++bound;
        } else {
            require(false, Code::Unsupported);
        }
    }
};

} // namespace

core::Result<void>
bind_observation_file(const SchemaCase& source, const ObservationFile& file,
                      const Exact& evaluation_time, structure::ExactLimits& exact,
                      std::uint64_t& work, core::ResourcePlan& plan, const std::stop_token& stop,
                      std::optional<ObservationBindingDiagnostic>& diagnostic) noexcept {
    diagnostic.reset();
    return core::guard_entry(operation, Code::InvalidInput, {}, [&]() -> core::Result<void> {
        ObservationBindingDiagnostic current{file.source, {}, {}, {}};
        try {
            require(!stop.stop_requested(), Code::Cancelled);
            const auto identity = take(source.identity());
            const auto& parent = file.source.parent;
            require(identity.source == parent.source &&
                        identity.contract_revision == parent.contract_revision &&
                        identity.representation_revision == parent.representation_revision,
                    Code::StaleGeneration);
            const auto root = take(source.root());
            const auto analysis = root.analysis.record();
            const auto* structural = std::get_if<schema::AnalysisStatic>(&analysis);
            require(structural && structural->physics.text() == "structural", Code::Unsupported);
            const auto units = root.units.record();
            Binding binding{file,
                            evaluation_time,
                            stop,
                            {static_cast<SiLengthUnit>(units.length.selected()),
                             static_cast<SiMassUnit>(units.mass.selected()),
                             static_cast<SiTimeUnit>(units.time.selected()),
                             static_cast<SiAngleUnit>(units.angle.selected())},
                            exact,
                            work,
                            plan,
                            current};
            try {
                for (core::Index i = 0; std::cmp_less(i, root.observations.size()); ++i) {
                    binding.diagnostic.observation.reset();
                    binding.diagnostic.case_index.reset();
                    binding.diagnostic.input_node.reset();
                    const auto observation = take(root.observations.at(i)).record();
                    std::visit([&](const auto& value) { binding.match(value, i); }, observation);
                }
                binding.step();
                require(binding.bound == file.records.size());
            } catch (const core::Failure&) {
                current = binding.diagnostic;
                throw;
            }
            exact = binding.exact;
            work = binding.work;
            return {};
        } catch (const core::Failure& failure) {
            diagnostic = current;
            return std::unexpected(failure);
        }
    });
}

} // namespace uvlas::io::detail
