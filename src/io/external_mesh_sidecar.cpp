#include "external_mesh_sidecar.hpp"

#include "json_object_reader.hpp"

#include <algorithm>
#include <charconv>
#include <initializer_list>
#include <limits>

namespace uvlas::io::detail {
namespace {

constexpr auto operation = core::Operation::PrepareModel;
using Code = core::FailureCode;
using Reason = MeshSidecarFailure;
using Q = structure::ExactRational;
using namespace std::string_view_literals;

struct Stop {
    core::Failure failure;
};

struct Reader {
    const MeshSidecarLimits& limits;
    core::ResourcePlan& plan;
    const std::stop_token& stop;
    RootSnapshotLocation source;
    std::optional<MeshSidecarDiagnostic>* diagnostic;
    std::optional<core::Index> node;
    std::uint64_t work{};
    std::uint64_t rows{};
    std::uint64_t faces{};

    [[noreturn]] void fail(Reason reason, Code code = Code::InvalidInput) const {
        if (diagnostic && code != Code::Cancelled) {
            *diagnostic = MeshSidecarDiagnostic{reason, source, node, {}, {}, {}, {}, {}};
        }
        throw Stop{{code, operation}};
    }

    void step(std::uint64_t amount = 1) {
        if (stop.stop_requested()) {
            fail(Reason::Resource, Code::Cancelled);
        }
        if (amount > limits.maximum_work - work) {
            fail(Reason::Resource, Code::ResourceExceeded);
        }
        work += amount;
    }

    void select(InputNodeView value) {
        node = value.ordinal().value();
        step();
    }

    template <class T>
    [[nodiscard]] T get(std::optional<T> value, Reason reason) const {
        if (!value) {
            fail(reason);
        }
        return *value;
    }

    template <class T>
    [[nodiscard]] T take(core::Result<T> value, Reason reason) const {
        if (!value) {
            if (diagnostic && value.error().code != Code::Cancelled) {
                *diagnostic = MeshSidecarDiagnostic{reason, source, node, {}, {}, {}, {}, {}};
            }
            throw Stop{value.error()};
        }
        return std::move(*value);
    }

    [[nodiscard]] InputObjectView object(InputNodeView value,
                                         std::initializer_list<std::string_view> required,
                                         std::initializer_list<std::string_view> optional = {}) {
        select(value);
        const auto obj = get(value.object(), Reason::Shape);
        for (core::Index i = 0; i < obj.size(); ++i) {
            const auto member = take(obj.at(i), Reason::Shape);
            select(member.value);
            step(member.key.size());
            if (std::ranges::find(required, member.key) == required.end() &&
                std::ranges::find(optional, member.key) == optional.end()) {
                fail(Reason::Key);
            }
        }
        select(value);
        for (const auto key : required) {
            step();
            if (!obj.find(key)) {
                fail(Reason::Key);
            }
        }
        return obj;
    }

    [[nodiscard]] InputNodeView field(const InputObjectView& object, std::string_view key) {
        const auto value = get(object.find(key), Reason::Key);
        select(value);
        return value;
    }

    [[nodiscard]] std::string_view string(InputNodeView value) {
        select(value);
        const auto text = get(value.string(), Reason::Scalar);
        step(text.size());
        return text;
    }

    [[nodiscard]] std::uint64_t id(InputNodeView value, bool positive = false) {
        const auto text = string(value);
        std::optional<ScalarConversionDiagnostic> detail;
        const auto result = read_id_digits(text, detail, {});
        if (!result || (positive && *result == 0)) {
            fail(Reason::Scalar);
        }
        return *result;
    }

    [[nodiscard]] std::uint64_t nastran_id(InputNodeView value) {
        const auto parsed = id(value, true);
        if (parsed > static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max())) {
            fail(Reason::Scalar);
        }
        return parsed;
    }

    [[nodiscard]] std::uint64_t integer(InputNodeView value, std::uint64_t maximum) {
        select(value);
        if (const auto u = value.unsigned_integer(); u && *u <= maximum) {
            return *u;
        }
        if (const auto i = value.signed_integer();
            i && *i >= 0 && static_cast<std::uint64_t>(*i) <= maximum) {
            return static_cast<std::uint64_t>(*i);
        }
        fail(Reason::Scalar);
    }

    template <class T, class Build>
    [[nodiscard]] core::PlannedBuffer<T> array(InputNodeView value, Build build,
                                               bool nonempty = false) {
        select(value);
        const auto array = get(value.array(), Reason::Shape);
        const auto count = static_cast<std::uint64_t>(array.size());
        if (nonempty && count == 0) {
            fail(Reason::Shape);
        }
        if (count > limits.maximum_rows - rows) {
            fail(Reason::Resource, Code::ResourceExceeded);
        }
        rows += count;
        auto result = take(plan.allocate<T>(array.size(), operation), Reason::Resource);
        for (core::Index i = 0; i < array.size(); ++i) {
            const auto item = take(array.at(i), Reason::Shape);
            select(item);
            result.values()[static_cast<std::size_t>(i)] = build(item);
        }
        return result;
    }

    [[nodiscard]] Q rational(InputNodeView value) {
        const auto row = object(value, {"numerator", "denominator"});
        const auto n = string(field(row, "numerator"));
        const auto d = string(field(row, "denominator"));
        if (d.starts_with('-')) {
            fail(Reason::Scalar);
        }
        return take(Q::from_integer_pair(n, d, limits.exact, plan), Reason::Scalar);
    }

    template <class Integer>
    [[nodiscard]] Q positive_integer(Integer value) {
        if (value <= 0) {
            fail(Reason::Scalar);
        }
        std::array<char, 24> text{};
        const auto digits = std::to_chars(text.data(), text.data() + text.size(), value);
        return take(Q::from_decimal(std::string_view{text.data(), digits.ptr}, limits.exact, plan),
                    Reason::Scalar);
    }

    [[nodiscard]] Q positive(InputNodeView value) {
        select(value);
        if (const auto decimal = value.decimal_token()) {
            step(decimal->size());
            auto exact = take(Q::from_decimal(*decimal, limits.exact, plan), Reason::Scalar);
            if (exact.sign() <= 0) {
                fail(Reason::Scalar);
            }
            return exact;
        }
        if (const auto integer = value.unsigned_integer()) {
            return positive_integer(*integer);
        }
        if (const auto integer = value.signed_integer()) {
            return positive_integer(*integer);
        }
        fail(Reason::Scalar);
    }

    [[nodiscard]] MeshSidecarChart chart(InputNodeView value) {
        const auto row = object(value, {"vertices", "orientation"});
        const auto word = string(field(row, "orientation"));
        if (word != "same" && word != "reversed") {
            fail(Reason::Chart);
        }
        const auto vertices_node = field(row, "vertices");
        const auto vertices = get(vertices_node.array(), Reason::Chart);
        if (vertices.size() != 3 && vertices.size() != 4) {
            fail(Reason::Chart);
        }
        if (faces == limits.maximum_faces) {
            fail(Reason::Resource, Code::ResourceExceeded);
        }
        ++faces;
        MeshSidecarChart output;
        output.count = static_cast<std::uint8_t>(vertices.size());
        output.orientation =
            word == "same" ? MeshSidecarOrientation::Same : MeshSidecarOrientation::Reversed;
        for (core::Index i = 0; i < vertices.size(); ++i) {
            const auto v = take(vertices.at(i), Reason::Chart);
            select(v);
            const auto coordinates = get(v.array(), Reason::Chart);
            if (coordinates.size() != 2) {
                fail(Reason::Chart);
            }
            output.vertices[static_cast<std::size_t>(i)] = {
                rational(take(coordinates.at(0), Reason::Chart)),
                rational(take(coordinates.at(1), Reason::Chart))};
        }
        return output;
    }
};

template <class Unit, std::size_t N>
[[nodiscard]] Unit unit(Reader& r, InputNodeView value,
                        const std::array<std::string_view, N>& words) {
    const auto text = r.string(value);
    const auto found = std::ranges::find(words, text);
    if (found == words.end()) {
        r.fail(Reason::Scalar);
    }
    return static_cast<Unit>(found - words.begin());
}

[[nodiscard]] SiLengthUnit length(Reader& r, InputNodeView value) {
    return unit<SiLengthUnit>(r, value, std::array{"m"sv, "mm"sv, "cm"sv, "km"sv, "in"sv, "ft"sv});
}

[[nodiscard]] SiUnits units(Reader& r, InputNodeView value) {
    const auto row = r.object(value, {"length", "mass", "time", "angle"});
    return {length(r, r.field(row, "length")),
            unit<SiMassUnit>(r, r.field(row, "mass"), std::array{"kg"sv, "g"sv, "t"sv, "lbm"sv}),
            unit<SiTimeUnit>(r, r.field(row, "time"), std::array{"s"sv, "ms"sv, "min"sv}),
            unit<SiAngleUnit>(r, r.field(row, "angle"), std::array{"rad"sv, "deg"sv})};
}

[[nodiscard]] GmshKey physical(Reader& r, InputNodeView value) {
    const auto row = r.object(value, {"dimension", "tag"});
    const auto dimension = r.integer(r.field(row, "dimension"), 2);
    const auto tag = r.integer(r.field(row, "tag"), std::numeric_limits<std::int32_t>::max());
    if (tag == 0) {
        r.fail(Reason::Scalar);
    }
    return {static_cast<std::uint8_t>(dimension), static_cast<std::int32_t>(tag)};
}

[[nodiscard]] GmshSidecarShell gmsh_shell(Reader& r, InputNodeView value) {
    const auto row =
        r.object(value, {"physical", "shell", "faces"}, {"material_layers", "full_thickness"});
    GmshSidecarShell out;
    out.node = value.ordinal().value();
    out.physical = physical(r, r.field(row, "physical"));
    if (out.physical.dimension != 2) {
        r.fail(Reason::Scalar);
    }
    out.shell = MeshShellId{r.id(r.field(row, "shell"))};
    out.faces = r.array<GmshSidecarFace>(
        r.field(row, "faces"),
        [&](InputNodeView v) {
            const auto face = r.object(v, {"element", "chart"});
            return GmshSidecarFace{r.id(r.field(face, "element"), true),
                                   r.chart(r.field(face, "chart")), v.ordinal().value()};
        },
        true);
    if (const auto layers = row.find("material_layers")) {
        out.material_layers = r.array<GmshSidecarLayer>(
            *layers,
            [&](InputNodeView v) {
                const auto layer = r.object(v, {"layer", "material"});
                return GmshSidecarLayer{core::LayerId{r.id(r.field(layer, "layer"))},
                                        MeshMaterialId{r.id(r.field(layer, "material"))},
                                        v.ordinal().value()};
            },
            true);
    }
    if (const auto thickness = row.find("full_thickness")) {
        out.full_thickness = r.positive(*thickness);
        out.full_thickness_node = thickness->ordinal().value();
    }
    return out;
}

[[nodiscard]] GmshSidecar gmsh(Reader& r, const InputObjectView& row) {
    GmshSidecar out;
    out.length = length(r, r.field(row, "length"));
    out.shells = r.array<GmshSidecarShell>(r.field(row, "shells"),
                                           [&](InputNodeView v) { return gmsh_shell(r, v); });
    out.supports = r.array<GmshSidecarSupport>(r.field(row, "supports"), [&](InputNodeView v) {
        const auto support = r.object(v, {"physical", "support"});
        return GmshSidecarSupport{physical(r, r.field(support, "physical")),
                                  MeshSupportId{r.id(r.field(support, "support"))},
                                  v.ordinal().value()};
    });
    out.joints = r.array<GmshSidecarJoint>(r.field(row, "joints"), [&](InputNodeView v) {
        const auto joint = r.object(v, {"first", "second", "joint"});
        return GmshSidecarJoint{physical(r, r.field(joint, "first")),
                                physical(r, r.field(joint, "second")),
                                MeshJointId{r.id(r.field(joint, "joint"))}, v.ordinal().value()};
    });
    return out;
}

template <SchemaNamespace Space>
[[nodiscard]] auto element_bindings(Reader& r, InputNodeView input, std::string_view name) {
    using Row = NastranSidecarElementBinding<Space>;
    return r.array<Row>(input, [&](InputNodeView v) {
        const auto row = r.object(v, {"element", name});
        return Row{NastranElementId{r.nastran_id(r.field(row, "element"))},
                   MeshCaseId<Space>{r.id(r.field(row, name))}, v.ordinal().value()};
    });
}

[[nodiscard]] std::uint8_t components(Reader& r, InputNodeView input) {
    const auto values = r.get(input.array(), Reason::Shape);
    if (values.size() < 1 || values.size() > 6) {
        r.fail(Reason::Scalar);
    }
    unsigned bits = 0;
    std::uint64_t last = 0;
    for (core::Index i = 0; i < values.size(); ++i) {
        const auto component = r.integer(r.take(values.at(i), Reason::Shape), 6);
        if (component <= last) {
            r.fail(Reason::Scalar);
        }
        bits |= 1U << (component - 1);
        last = component;
    }
    return static_cast<std::uint8_t>(bits);
}

[[nodiscard]] NastranSidecar nastran(Reader& r, const InputObjectView& row) {
    NastranSidecar out;
    const auto generator = r.string(r.field(row, "generator"));
    if (generator.empty() || !std::in_range<core::Index>(generator.size())) {
        r.fail(Reason::Scalar);
    }
    out.generator =
        r.take(r.plan.allocate<char>(static_cast<core::Index>(generator.size()), operation),
               Reason::Resource);
    std::ranges::copy(generator, out.generator.values().begin());
    out.units = units(r, r.field(row, "units"));
    out.shells = r.array<NastranSidecarShell>(r.field(row, "shells"), [&](InputNodeView v) {
        const auto shell = r.object(v, {"element", "shell", "chart"});
        return NastranSidecarShell{NastranElementId{r.nastran_id(r.field(shell, "element"))},
                                   MeshShellId{r.id(r.field(shell, "shell"))},
                                   r.chart(r.field(shell, "chart")), v.ordinal().value()};
    });
    out.materials =
        r.array<NastranSidecarMaterial>(r.field(row, "materials"), [&](InputNodeView v) {
            const auto material = r.object(v, {"mid", "material"});
            return NastranSidecarMaterial{NastranMaterialId{r.nastran_id(r.field(material, "mid"))},
                                          MeshMaterialId{r.id(r.field(material, "material"))},
                                          v.ordinal().value()};
        });
    out.beams = element_bindings<SchemaNamespace::Beam>(r, r.field(row, "beams"), "beam");
    out.masses = element_bindings<SchemaNamespace::Mass>(r, r.field(row, "masses"), "mass");
    out.supports = r.array<NastranSidecarSupport>(r.field(row, "supports"), [&](InputNodeView v) {
        const auto support = r.object(v, {"set", "grid", "components", "support"});
        return NastranSidecarSupport{NastranSupportSetId{r.nastran_id(r.field(support, "set"))},
                                     NastranGridId{r.nastran_id(r.field(support, "grid"))},
                                     components(r, r.field(support, "components")),
                                     MeshSupportId{r.id(r.field(support, "support"))},
                                     v.ordinal().value()};
    });
    out.joints = element_bindings<SchemaNamespace::Joint>(r, r.field(row, "joints"), "joint");
    return out;
}

[[nodiscard]] std::optional<MeshConversionBudget> conversion(Reader& r,
                                                             const InputObjectView& root) {
    const auto value = root.find("conversion");
    if (!value) {
        return {};
    }
    const auto row =
        r.object(*value, {"distance_budget", "normal_budget", "conversion_observations"});
    return MeshConversionBudget{
        r.positive(r.field(row, "distance_budget")), r.positive(r.field(row, "normal_budget")),
        r.array<MeshConversionObservation>(
            r.field(row, "conversion_observations"),
            [&](InputNodeView v) {
                return MeshConversionObservation{core::ObservationId{r.id(v)}, v.ordinal().value()};
            },
            true)};
}

} // namespace

core::Result<ParsedMeshSidecar>
read_mesh_sidecar(const AcquiredRootSnapshot& snapshot, const MeshSidecarLimits& limits,
                  core::ResourcePlan& plan, const std::stop_token& stop,
                  std::optional<MeshSidecarDiagnostic>* diagnostic) noexcept {
    if (diagnostic) {
        diagnostic->reset();
    }
    return core::guard_entry(
        operation, Code::InvalidInput, {}, [&]() -> core::Result<ParsedMeshSidecar> {
            try {
                Reader reader{limits, plan, stop, snapshot.location, diagnostic, {}};
                reader.step();
                const auto file = reader.get(snapshot.file.view(), Reason::Snapshot);
                const bool is_gmsh = snapshot.location.role == ReferenceRole::GmshSidecar;
                if (!is_gmsh && snapshot.location.role != ReferenceRole::NastranSidecar) {
                    reader.fail(Reason::Snapshot);
                }
                std::optional<InputDiagnostic> json;
                auto doc =
                    read_json_object(file.bytes, snapshot.location.source, limits.json, plan, json);
                reader.step();
                if (!doc) {
                    if (diagnostic) {
                        *diagnostic = MeshSidecarDiagnostic{
                            Reason::Json, snapshot.location,
                            {},           json ? json->byte_offset : std::nullopt,
                            {},           {},
                            {},           {}};
                    }
                    auto error = doc.error();
                    if (error.code == Code::Nonfinite) {
                        error.code = Code::InvalidInput;
                    }
                    return std::unexpected(error);
                }
                const auto view = reader.take(doc->view(), Reason::Json);
                if (view.identity.digest != file.digest) {
                    reader.fail(Reason::Snapshot);
                }
                const auto envelope = reader.get(view.root.object(), Reason::Shape);
                const auto version = reader.string(reader.field(envelope, "schema"));
                if (version != (is_gmsh ? "uvlas.gmsh-sidecar/1" : "uvlas.nastran-sidecar/1")) {
                    const auto family = is_gmsh ? "uvlas.gmsh-sidecar/" : "uvlas.nastran-sidecar/";
                    reader.fail(Reason::Schema, version.starts_with(family) ? Code::Unsupported
                                                                            : Code::InvalidInput);
                }
                const auto body =
                    is_gmsh ? reader.object(
                                  view.root,
                                  {"schema", "geometry", "length", "shells", "supports", "joints"},
                                  {"conversion"})
                            : reader.object(view.root,
                                            {"schema", "geometry", "generator", "units", "shells",
                                             "materials", "beams", "masses", "supports", "joints"},
                                            {"conversion"});
                const MeshGeometryId geometry{reader.id(reader.field(body, "geometry"))};
                auto payload = is_gmsh ? MeshSidecarPayload{gmsh(reader, body)}
                                       : MeshSidecarPayload{nastran(reader, body)};
                auto budget = conversion(reader, body);
                reader.step();
                return ParsedMeshSidecar{snapshot.location,  file.digest,       geometry,
                                         std::move(payload), std::move(budget), reader.work};
            } catch (const Stop& failed) {
                return std::unexpected(failed.failure);
            }
        });
}
} // namespace uvlas::io::detail
