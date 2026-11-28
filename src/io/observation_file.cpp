#include "observation_file.hpp"

#include "json_object_reader.hpp"
#include "source_digest.hpp"

#include <algorithm>
#include <charconv>
#include <limits>

namespace uvlas::io::detail {
namespace {

using Code = core::FailureCode;
using Exact = structure::ExactRational;
constexpr auto operation = core::Operation::PrepareModel;

struct Reader {
    ObservationFileDiagnostic where;
    const std::stop_token& stop;
    ObservationFileLimits limits;
    core::ResourcePlan& plan;

    template <class T>
    T take(core::Result<T> value) {
        if (!value)
            throw value.error();
        return std::move(*value);
    }

    void require(bool condition, Code code = Code::InvalidInput) {
        if (!condition)
            throw core::Failure{code, operation};
    }

    void step(std::uint64_t work) {
        require(!stop.stop_requested(), Code::Cancelled);
        require(work <= limits.maximum_work, Code::ResourceExceeded);
        limits.maximum_work -= work;
    }

    void at(InputNodeView node) {
        where.node = node.ordinal();
        step(1);
    }

    InputObjectView object(InputNodeView node) {
        at(node);
        const auto object = node.object();
        if (!object)
            throw core::Failure{Code::InvalidInput, operation};
        return *object;
    }

    InputObjectView object(InputNodeView node, core::Index fields) {
        const auto value = object(node);
        require(value.size() == fields);
        return value;
    }

    InputArrayView array(InputNodeView node) {
        at(node);
        const auto array = node.array();
        if (!array)
            throw core::Failure{Code::InvalidInput, operation};
        return *array;
    }

    InputNodeView field(const InputObjectView& object, std::string_view key) {
        const auto fields = static_cast<std::uint64_t>(object.size());
        require(fields <= std::numeric_limits<std::uint64_t>::max() / (1 + key.size()),
                Code::ResourceExceeded);
        step(fields * (1 + key.size()));
        const auto node = object.find(key);
        if (!node)
            throw core::Failure{Code::InvalidInput, operation};
        at(*node);
        return *node;
    }

    std::string_view string(InputNodeView node) {
        at(node);
        const auto text = node.string();
        if (!text)
            throw core::Failure{Code::InvalidInput, operation};
        step(text->size());
        return *text;
    }

    std::string_view field_string(const InputObjectView& object, std::string_view key) {
        return string(field(object, key));
    }

    std::uint64_t id(const InputObjectView& object, std::string_view key) {
        const auto text = field_string(object, key);
        std::optional<ScalarConversionDiagnostic> diagnostic;
        return take(read_id_digits(text, diagnostic, {}));
    }

    Exact number(InputNodeView node) {
        at(node);
        step(limits.exact.maximum_work_bits);
        auto value = [&]() -> Exact {
            if (const auto decimal = node.decimal_token()) {
                step(decimal->size());
                return take(Exact::from_decimal(*decimal, limits.exact, plan));
            }
            const auto integer = node.native_integer();
            if (!integer)
                throw core::Failure{Code::InvalidInput, operation};
            return std::visit(
                [&](auto number) {
                    std::array<char, 24> buffer{};
                    const auto end =
                        std::to_chars(buffer.data(), buffer.data() + buffer.size(), number);
                    step(static_cast<std::size_t>(end.ptr - buffer.data()));
                    return take(Exact::from_integer_pair(
                        {buffer.data(), static_cast<std::size_t>(end.ptr - buffer.data())}, "1",
                        limits.exact, plan));
                },
                integer->selection());
        }();
        // Cumulative admission overestimates retained values by including each operation's
        // temporary limb forecast. Each following operation sees only the remaining allowance.
        if (auto* memory =
                std::get_if<structure::ExactEstimateAdmission>(&limits.exact.external_memory)) {
            const auto used = value.preparation().estimated_external_payload_bytes;
            require(used <= memory->bytes, Code::ResourceExceeded);
            memory->bytes -= used;
        }
        return value;
    }

    std::array<Exact, 3> vector(InputNodeView node) {
        const auto values = array(node);
        require(values.size() == 3);
        std::array<Exact, 3> output;
        for (core::Index i = 0; i != 3; ++i)
            output[static_cast<std::size_t>(i)] = number(take(values.at(i)));
        return output;
    }

    void world(const InputObjectView& descriptor) {
        const auto basis = object(field(descriptor, "basis"), 1);
        require(field_string(basis, "kind") == "world");
    }

    ObservationFileRecord record(InputNodeView input) {
        const auto value = object(input, 4);
        ObservationFileRecord result;
        result.input_node = input.ordinal();
        result.id = core::ObservationId{id(value, "id")};
        const auto descriptor_node = field(value, "observation");
        const auto descriptor = object(descriptor_node, 4);
        const auto kind = field_string(descriptor, "kind");
        world(descriptor);
        const bool point = kind == "point";
        if (point) {
            require(field_string(descriptor, "quantity") == "displacement");
            const auto physical = object(field(descriptor, "point"), 3);
            const auto target = object(field(physical, "target"), 2);
            require(field_string(target, "kind") == "shell");
            const auto shell = id(target, "id");
            const auto side = field_string(physical, "side");
            require(side == "interior" || side == "lower" || side == "upper");
            const auto location = side == "lower"   ? ObservationPointSide::Lower
                                  : side == "upper" ? ObservationPointSide::Upper
                                                    : ObservationPointSide::Interior;
            result.descriptor =
                ObservationFilePoint{shell, location, vector(field(physical, "coordinates"))};
        } else {
            require(kind == "reaction");
            const auto support = id(descriptor, "support");
            result.descriptor =
                ObservationFileReaction{support, vector(field(descriptor, "origin_m"))};
        }
        const auto comparison = object(field(value, "comparison"), 3);
        require(field_string(comparison, "statistic") == "instantaneous");
        require(field_string(comparison, "time_side") == "right");
        result.time_seconds = number(field(comparison, "time_s"));

        const auto components = array(field(value, "components"));
        constexpr std::array<std::string_view, 3> xyz{"x", "y", "z"};
        constexpr std::array<std::string_view, 6> wrench{"Fx", "Fy", "Fz", "Mx", "My", "Mz"};
        const auto names = point ? std::span<const std::string_view>{xyz}
                                 : std::span<const std::string_view>{wrench};
        require(components.size() > 0 && std::cmp_less_equal(components.size(), names.size()));
        for (core::Index i = 0; i < components.size(); ++i) {
            const auto element = object(take(components.at(i)), 4);
            const auto name = field_string(element, "name");
            const auto index =
                static_cast<std::size_t>(std::ranges::find(names, name) - names.begin());
            require(index < names.size() && !result.components[index]);
            const std::string_view unit = point ? "m" : index < 3 ? "N" : "N*m";
            require(field_string(element, "unit") == unit);
            auto center = number(field(element, "value"));
            auto uncertainty = number(field(element, "uncertainty"));
            require(uncertainty.sign() >= 0);
            result.components[index].emplace(
                ObservationFileScalar{std::move(center), std::move(uncertainty)});
        }
        return result;
    }
};

bool same_parent(const SchemaIdentity& a, const SchemaIdentity& b) {
    return a.source == b.source && a.contract_revision == b.contract_revision &&
           a.representation_revision == b.representation_revision;
}

} // namespace

core::Result<ObservationFile>
read_observation_file(const AcquiredRootSnapshot& snapshot, const SchemaIdentity& parent,
                      const SourceDigest& condition, std::span<const std::uint64_t> expected,
                      ObservationFileLimits& limits, core::ResourcePlan& plan,
                      const std::stop_token& stop,
                      std::optional<ObservationFileDiagnostic>& diagnostic) noexcept {
    diagnostic.reset();
    return core::guard_entry(
        operation, Code::InvalidInput, {}, [&]() -> core::Result<ObservationFile> {
            Reader read{{snapshot.location, {}, {}}, stop, limits, plan};
            try {
                read.step(1);
                read.require(same_parent(snapshot.location.parent, parent), Code::StaleGeneration);
                read.require(snapshot.location.role == ReferenceRole::ObservationReference);
                const auto file = snapshot.file.view();
                read.require(file.has_value());
                std::optional<InputDiagnostic> syntax;
                read.step(file->bytes.size());
                const auto document = read_json_object(file->bytes, snapshot.location.source,
                                                       limits.json, plan, syntax);
                read.step(1);
                if (!document) {
                    if (syntax)
                        read.where.byte_offset = syntax->byte_offset;
                    auto error = document.error();
                    if (error.code == Code::Nonfinite)
                        error.code = Code::InvalidInput;
                    throw error;
                }
                const auto input = read.take(document->view());
                read.require(input.identity.digest == file->digest);
                const auto object = read.object(input.root);
                read.require(read.field_string(object, "schema") == "uvlas.observation-reference/1",
                             Code::Unsupported);
                read.at(input.root);
                read.require(object.size() == 3);
                const auto digest =
                    decode_source_digest(read.field_string(object, "condition_sha256"));
                read.require(digest && *digest == condition);
                const auto entries = read.array(read.field(object, "observations"));
                read.require(entries.size() > 0 && std::cmp_equal(entries.size(), expected.size()));
                auto output =
                    read.take(plan.allocate<ObservationFileRecord>(entries.size(), operation));
                for (core::Index i = 0; i < entries.size(); ++i) {
                    read.step(static_cast<std::uint64_t>(i) + 1);
                    auto record = read.record(read.take(entries.at(i)));
                    read.step(expected.size());
                    read.require(std::ranges::find(expected, record.id.value()) != expected.end());
                    output.values()[static_cast<std::size_t>(i)] = std::move(record);
                }
                std::ranges::sort(output.values(), [&](const auto& a, const auto& b) {
                    read.step(1);
                    return a.id.value() < b.id.value();
                });
                const auto pair =
                    std::ranges::adjacent_find(output.values(), [&](const auto& a, const auto& b) {
                        read.step(1);
                        return a.id == b.id;
                    });
                if (pair != output.values().end()) {
                    read.where.node = pair->input_node;
                    read.require(false);
                }
                read.step(1);
                limits = read.limits;
                return ObservationFile{snapshot.location, file->digest, *digest, std::move(output)};
            } catch (const core::Failure& failure) {
                diagnostic = read.where;
                return std::unexpected(failure);
            }
        });
}

} // namespace uvlas::io::detail
