#include "visualization_series.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <limits>
#include <tuple>
#include <type_traits>
#include <utility>

namespace uvlas::io::detail {
namespace {
constexpr auto operation = core::Operation::WriteResults;
constexpr auto invalid = core::FailureCode::InvalidInput;
constexpr auto capacity = core::FailureCode::ResourceExceeded;
constexpr auto nonfinite = core::FailureCode::Nonfinite;

void require(bool condition, core::FailureCode code = invalid) {
    if (!condition)
        throw core::Failure{code, operation};
}

void poll(const std::stop_token& stop) {
    require(!stop.stop_requested(), core::FailureCode::Cancelled);
}

void count(std::size_t& total, std::size_t add, std::size_t maximum) {
    require(add <= maximum - total, capacity);
    total += add;
}

template <class T>
T take(core::Result<T> result) {
    if (!result)
        throw result.error();
    return std::move(*result);
}

template <class T>
auto allocate(core::ResourcePlan& plan, std::size_t size) {
    require(std::in_range<core::Index>(size), capacity);
    return take(plan.allocate<T>(static_cast<core::Index>(size), operation));
}

template <class Scalar>
struct NumberText {
    explicit NumberText(Scalar value) {
        const auto [end, error] = std::to_chars(bytes.data(), bytes.data() + bytes.size(), value);
        require(error == std::errc{}, capacity);
        size = static_cast<std::size_t>(end - bytes.data());
    }
    [[nodiscard]] std::string_view view() const noexcept {
        return {bytes.data(), size};
    }

    std::array<char, 32> bytes{};
    std::size_t size = 0;
};

struct FileName {
    void append(std::string_view text) {
        require(text.size() <= bytes.size() - size, capacity);
        std::ranges::copy(text, std::span{bytes}.subspan(size).begin());
        size += text.size();
    }
    [[nodiscard]] std::string_view view() const noexcept {
        return {bytes.data(), size};
    }

    std::array<char, 64> bytes{};
    std::size_t size = 0;
};

FileName vtu_name(core::AcceptedStateId id, std::size_t part) {
    FileName name;
    name.append("accepted-");
    name.append(NumberText{id.value()}.view());
    name.append("-part-");
    name.append(NumberText{part}.view());
    name.append(".vtu");
    return name;
}

FileName pvd_name(std::uint64_t segment) {
    FileName name;
    name.append("segment-");
    name.append(NumberText{segment}.view());
    name.append(".pvd");
    return name;
}

class Xml {
public:
    Xml(std::size_t maximum_bytes, const std::stop_token& stop,
        std::span<char> storage = {}) noexcept
        : maximum_(maximum_bytes), storage_(storage), stop_(stop) {}

    void append(std::string_view text) {
        require(text.size() <= maximum_ - size_, capacity);
        if (!storage_.empty())
            std::ranges::copy(text, storage_.subspan(size_).begin());
        size_ += text.size();
    }

    template <class Scalar>
    void number(Scalar value) {
        append(NumberText{value}.view());
    }

    template <class Scalar>
    void item(Scalar value) {
        number(value);
        append(" ");
    }

    template <class Vector>
    void tuple(const Vector& values) {
        for (const auto value : values)
            item(value);
    }

    [[nodiscard]] std::size_t size() const noexcept {
        return size_;
    }
    void point() const {
        poll(stop_);
    }

private:
    std::size_t maximum_;
    std::span<char> storage_;
    const std::stop_token& stop_;
    std::size_t size_ = 0;
};

constexpr std::array<std::string_view, 4> measures{"green", "second_pk", "first_pk", "cauchy"};
constexpr std::array<std::string_view, 2> bases{"world", "material"};
constexpr std::array<std::string_view, 3> contributions{"total", "elastic", "viscous"};

template <class Vector>
bool finite(const Vector& v) noexcept {
    return std::ranges::all_of(v, [](double value) { return std::isfinite(value); });
}

bool nonzero(const std::array<double, 3>& vector) noexcept {
    return std::ranges::any_of(vector, [](double value) { return value != 0.0; });
}

void validate_point(const VisualizationPoint& value, VisualizationLayerSide side) {
    require(finite(value.reference_metres) && finite(value.current_metres) &&
                finite(value.reference_normal) && finite(value.current_normal) &&
                std::isfinite(value.reference_thickness_metres) &&
                std::isfinite(value.current_thickness_metres) &&
                std::isfinite(value.layer_fraction),
            nonfinite);
    require(value.reference_thickness_metres > 0.0 && value.current_thickness_metres > 0.0 &&
            nonzero(value.reference_normal) && nonzero(value.current_normal));
    switch (side) {
    case VisualizationLayerSide::Interior:
        require(value.layer_fraction > 0.0 && value.layer_fraction < 1.0);
        return;
    case VisualizationLayerSide::Lower:
        require(value.layer_fraction == 0.0);
        return;
    case VisualizationLayerSide::Upper:
        require(value.layer_fraction == 1.0);
        return;
    }
    require(false);
}

void validate(const VisualizationPart& part, const std::stop_token& stop) {
    require(!part.points.empty() && !part.triangles.empty());
    constexpr auto index_max = static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max());
    require(part.points.size() <= index_max && part.triangles.size() <= index_max / 3, capacity);
    for (const auto& point : part.points) {
        poll(stop);
        validate_point(point, part.key.side);
    }
    for (const auto& triangle : part.triangles) {
        poll(stop);
        require(triangle[0] < part.points.size() && triangle[1] < part.points.size() &&
                triangle[2] < part.points.size() && triangle[0] != triangle[1] &&
                triangle[0] != triangle[2] && triangle[1] != triangle[2]);
    }

    using Contributions = std::array<bool, contributions.size()>;
    using Bases = std::array<Contributions, bases.size()>;
    std::array<Bases, measures.size()> seen{};
    require(part.tensors.size() <= measures.size() * bases.size() * contributions.size());
    for (const auto& tensor : part.tensors) {
        const auto measure = static_cast<std::size_t>(tensor.measure);
        const auto basis = static_cast<std::size_t>(tensor.basis);
        const auto contribution = static_cast<std::size_t>(tensor.contribution);
        require(measure < measures.size() && basis < bases.size() &&
                contribution < contributions.size());
        require(tensor.measure != VisualizationMeasure::Green ||
                tensor.contribution == VisualizationContribution::Total);
        auto& exists = seen[measure][basis][contribution];
        require(!exists);
        exists = true;
        require(std::holds_alternative<std::span<const VisualizationMatrix>>(tensor.values) ==
                (tensor.measure == VisualizationMeasure::FirstPk));
        std::visit(
            [&](auto values) {
                require(values.size() == part.points.size());
                for (const auto& value : values) {
                    poll(stop);
                    require(finite(value), nonfinite);
                }
            },
            tensor.values);
    }
}

std::size_t validate(std::span<const VisualizationFrame> frames, VisualizationLimits limits,
                     const std::stop_token& stop) {
    require(frames.size() <= limits.maximum_frames, capacity);
    std::size_t files = 0;
    std::size_t points = 0;
    std::size_t triangles = 0;
    const VisualizationFrame* previous = nullptr;
    for (const auto& frame : frames) {
        poll(stop);
        require(std::isfinite(frame.time_seconds), nonfinite);
        require(!frame.parts.empty());
        require(!previous ||
                (previous->accepted < frame.accepted && previous->segment <= frame.segment));
        const bool same_segment = previous && previous->segment == frame.segment;
        if (previous && same_segment) {
            require(previous->time_seconds < frame.time_seconds &&
                    previous->parts.size() == frame.parts.size());
        } else {
            count(files, 1, limits.maximum_files);
        }
        count(files, frame.parts.size(), limits.maximum_files);
        for (std::size_t i = 0; i < frame.parts.size(); ++i) {
            const auto& part = frame.parts[i];
            require(i == 0 || frame.parts[i - 1].key < part.key);
            if (previous && same_segment)
                require(previous->parts[i].key == part.key);
            count(points, part.points.size(), limits.maximum_points);
            count(triangles, part.triangles.size(), limits.maximum_triangles);
            validate(part, stop);
        }
        previous = &frame;
    }
    return files;
}

constexpr std::string_view xml_start = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
constexpr std::string_view array_end = "</DataArray>\n";

void open_array(Xml& xml, std::string_view type, std::string_view name) {
    xml.append("<DataArray type=\"");
    xml.append(type);
    xml.append("\" Name=\"");
    xml.append(name);
}

template <class Scalar>
void field(Xml& xml, std::string_view type, std::string_view name, Scalar value) {
    open_array(xml, type, name);
    xml.append("\" NumberOfTuples=\"1\" format=\"ascii\">");
    xml.number(value);
    xml.append(array_end);
}

template <class Member>
void point_data(Xml& xml, const VisualizationPart& part, std::string_view name, Member member) {
    open_array(xml, "Float64", name);
    if constexpr (!std::is_same_v<Member, double VisualizationPoint::*>)
        xml.append("\" NumberOfComponents=\"3");
    xml.append("\" format=\"ascii\">");
    for (const auto& point : part.points) {
        xml.point();
        if constexpr (std::is_same_v<Member, double VisualizationPoint::*>)
            xml.item(point.*member);
        else
            xml.tuple(point.*member);
    }
    xml.append(array_end);
}

void tensor_data(Xml& xml, const VisualizationTensorField& tensor) {
    const bool green = tensor.measure == VisualizationMeasure::Green;
    open_array(xml, "Float64", green ? "strain_" : "stress_");
    xml.append(measures[static_cast<std::size_t>(tensor.measure)]);
    if (!green) {
        xml.append("_");
        xml.append(contributions[static_cast<std::size_t>(tensor.contribution)]);
    }
    xml.append("_");
    xml.append(bases[static_cast<std::size_t>(tensor.basis)]);
    xml.append(green ? "_1" : "_Pa");
    xml.append("\" NumberOfComponents=\"9\" format=\"ascii\">");
    std::visit(
        [&](auto values) {
            for (const auto& value : values) {
                xml.point();
                if constexpr (std::tuple_size_v<std::remove_cvref_t<decltype(value)>> == 6) {
                    constexpr std::array<std::size_t, 9> order{0, 5, 4, 5, 1, 3, 4, 3, 2};
                    for (const auto index : order)
                        xml.item(value[index]);
                } else {
                    xml.tuple(value);
                }
            }
        },
        tensor.values);
    xml.append(array_end);
}

void vtu(Xml& xml, const VisualizationFrame& frame, const VisualizationPart& part) {
    xml.append(xml_start);
    xml.append("<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n"
               "<UnstructuredGrid>\n<FieldData>\n");
    field(xml, "Float64", "time_s", frame.time_seconds);
    field(xml, "UInt64", "accepted_state_id", frame.accepted.value());
    field(xml, "UInt64", "segment", frame.segment);
    field(xml, "UInt64", "component_id", part.key.component.value());
    field(xml, "UInt64", "layer_id", part.key.layer.value());
    field(xml, "UInt8", "layer_side", static_cast<unsigned>(part.key.side));
    xml.append("</FieldData>\n<Piece NumberOfPoints=\"");
    xml.number(part.points.size());
    xml.append("\" NumberOfCells=\"");
    xml.number(part.triangles.size());
    xml.append("\">\n<PointData Normals=\"current_normal\">\n");
    point_data(xml, part, "reference_position_m", &VisualizationPoint::reference_metres);
    point_data(xml, part, "reference_normal", &VisualizationPoint::reference_normal);
    point_data(xml, part, "current_normal", &VisualizationPoint::current_normal);
    point_data(xml, part, "reference_thickness_m", &VisualizationPoint::reference_thickness_metres);
    point_data(xml, part, "current_thickness_m", &VisualizationPoint::current_thickness_metres);
    point_data(xml, part, "layer_fraction", &VisualizationPoint::layer_fraction);
    for (const auto& tensor : part.tensors)
        tensor_data(xml, tensor);
    xml.append("</PointData>\n<CellData/>\n<Points>\n");
    point_data(xml, part, "current_position_m", &VisualizationPoint::current_metres);
    xml.append("</Points>\n<Cells>\n"
               "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">");
    for (const auto& triangle : part.triangles) {
        xml.point();
        xml.tuple(triangle);
    }
    xml.append(array_end);
    xml.append("<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">");
    for (std::size_t i = 0; i < part.triangles.size(); ++i) {
        xml.point();
        xml.item((i + 1) * 3);
    }
    xml.append(array_end);
    xml.append("<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">");
    for (std::size_t i = 0; i < part.triangles.size(); ++i) {
        xml.point();
        xml.item(5); // VTK_TRIANGLE
    }
    xml.append(array_end);
    xml.append("</Cells>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n");
}

void pvd(Xml& xml, std::span<const VisualizationFrame> segment) {
    xml.append(xml_start);
    xml.append("<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\">\n"
               "<Collection>\n");
    for (const auto& frame : segment) {
        for (std::size_t part = 0; part < frame.parts.size(); ++part) {
            xml.point();
            xml.append("<DataSet timestep=\"");
            xml.number(frame.time_seconds);
            xml.append("\" group=\"\" part=\"");
            xml.number(part);
            xml.append("\" file=\"");
            xml.append(vtu_name(frame.accepted, part).view());
            xml.append("\"/>\n");
        }
    }
    xml.append("</Collection>\n</VTKFile>\n");
}

template <class Record>
void write_series(Xml& xml, std::span<const VisualizationFrame> frames, Record record) {
    while (!frames.empty()) {
        const auto id = frames.front().segment;
        std::size_t length = 1;
        while (length < frames.size() && frames[length].segment == id) {
            xml.point();
            ++length;
        }
        const auto segment = frames.first(length);
        auto begin = xml.size();
        pvd(xml, segment);
        record(pvd_name(id), begin, xml.size() - begin);
        for (const auto& frame : segment) {
            for (std::size_t part = 0; part < frame.parts.size(); ++part) {
                begin = xml.size();
                vtu(xml, frame, frame.parts[part]);
                record(vtu_name(frame.accepted, part), begin, xml.size() - begin);
            }
        }
        frames = frames.subspan(segment.size());
    }
}

} // namespace

VisualizationSeries::VisualizationSeries(core::PlannedBuffer<FileLocation> files,
                                         core::PlannedBuffer<char> bytes) noexcept
    : files_(std::move(files)), bytes_(std::move(bytes)) {}

std::size_t VisualizationSeries::size() const noexcept {
    return files_.size();
}

std::optional<VisualizationFileView> VisualizationSeries::at(std::size_t index) const& noexcept {
    if (index >= files_.size())
        return std::nullopt;
    const auto& file = files_.values()[index];
    const auto text = bytes_.values().subspan(file.offset, file.size);
    return VisualizationFileView{{file.name.data(), file.name_size}, {text.data(), text.size()}};
}

core::Result<VisualizationSeries>
encode_visualization_series(std::span<const VisualizationFrame> frames, VisualizationLimits limits,
                            core::ResourcePlan& plan, const std::stop_token& stop) noexcept {
    return core::guard_entry(operation, invalid, {}, [&]() -> core::Result<VisualizationSeries> {
        try {
            poll(stop);
            const auto n = validate(frames, limits, stop);
            using Location = VisualizationSeries::FileLocation;
            require(std::in_range<core::Index>(n), capacity);
            const auto metadata = take(
                core::checked_byte_count(static_cast<core::Index>(n), sizeof(Location), operation));
            require(metadata <= plan.remaining_bytes(), capacity);
            const auto payload = std::min(static_cast<std::uint64_t>(limits.maximum_output_bytes),
                                          plan.remaining_bytes() - metadata);
            Xml measure{static_cast<std::size_t>(payload), stop};
            write_series(measure, frames, [](FileName, std::size_t, std::size_t) {});
            require(std::in_range<core::Index>(measure.size()), capacity);
            poll(stop);
            auto files = allocate<Location>(plan, n);
            auto bytes = allocate<char>(plan, measure.size());
            Xml output{bytes.size(), stop, bytes.values()};
            std::size_t index = 0;
            write_series(output, frames, [&](FileName name, std::size_t offset, std::size_t size) {
                files.values()[index++] = {name.bytes, name.size, offset, size};
            });
            poll(stop);
            return VisualizationSeries{std::move(files), std::move(bytes)};
        } catch (const core::Failure& error) {
            return std::unexpected(error);
        }
    });
}

} // namespace uvlas::io::detail
