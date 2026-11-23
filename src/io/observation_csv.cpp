#include "observation_csv.hpp"

#include "utf8.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>

namespace uvlas::io::detail {
namespace {

constexpr auto operation = core::Operation::WriteResults;
constexpr auto invalid = core::FailureCode::InvalidInput;
constexpr auto capacity = core::FailureCode::ResourceExceeded;
constexpr auto nonfinite = core::FailureCode::Nonfinite;
constexpr std::string_view header =
    "time_hi_s,time_lo_s,checkpoint,segment,observation,component,quantity,unit,basis,"
    "reference_point,layer,side,representation,value_hi,value_lo,radius,exponent2\r\n";

void require(bool condition, core::FailureCode code = invalid) {
    if (!condition)
        throw core::Failure{code, operation};
}

void poll(const std::stop_token& stop) {
    require(!stop.stop_requested(), core::FailureCode::Cancelled);
}

template <class T>
[[nodiscard]] auto allocate(core::ResourcePlan& resources, std::size_t size) {
    require(std::in_range<core::Index>(size), capacity);
    auto storage = resources.allocate<T>(static_cast<core::Index>(size), operation);
    if (!storage)
        throw storage.error();
    return std::move(*storage);
}

[[nodiscard]] std::string_view side_name(ObservationSide side) noexcept {
    switch (side) {
    case ObservationSide::None:
        return "none";
    case ObservationSide::Interior:
        return "interior";
    case ObservationSide::Lower:
        return "lower";
    case ObservationSide::Upper:
        return "upper";
    }
    return {};
}

struct ColumnKey {
    core::ObservationId observation{0};
    std::uint64_t component = 0;
    std::optional<core::LayerId> layer;
    ObservationSide side = ObservationSide::None;

    auto operator<=>(const ColumnKey&) const = default;
};

void validate_columns(std::span<const ObservationCsvColumn> columns, ObservationCsvLimits limits,
                      core::ResourcePlan& resources, const std::stop_token& stop) {
    auto keys = allocate<ColumnKey>(resources, columns.size());
    for (std::size_t i = 0; i < columns.size(); ++i) {
        poll(stop);
        const auto& column = columns[i];
        require(!side_name(column.side).empty());
        require(column.layer.has_value() == (column.side != ObservationSide::None));
        for (const auto label :
             {column.quantity, column.unit, column.basis, column.reference_point}) {
            require(label.size() <= limits.maximum_label_bytes, capacity);
            require(!label.empty() && label.find('\0') == label.npos && valid_utf8(label));
        }
        keys.values()[i] = {column.observation, column.component, column.layer, column.side};
    }
    std::ranges::sort(keys.values());
    require(std::ranges::adjacent_find(keys.values()) == keys.values().end());
}

class CsvWriter {
public:
    explicit CsvWriter(std::size_t limit) noexcept : limit_(limit) {}
    explicit CsvWriter(std::span<char> output) noexcept : limit_(output.size()), output_(output) {}

    void append(std::string_view text) {
        require(text.size() <= limit_ - size_, capacity);
        if (!output_.empty())
            std::ranges::copy(text, output_.subspan(size_, text.size()).begin());
        size_ += text.size();
    }

    void label(std::string_view text) {
        if (text.find_first_of(",\"\r\n") == text.npos) {
            append(text);
            return;
        }
        append("\"");
        for (auto quote = text.find('"'); quote != text.npos; quote = text.find('"')) {
            append(text.substr(0, quote));
            append("\"\"");
            text.remove_prefix(quote + 1);
        }
        append(text);
        append("\"");
    }

    template <class Number>
    void number(Number value) {
        // Covers signed/unsigned 64-bit integers and shortest-round-trip binary64.
        std::array<char, 32> buffer{};
        const auto [end, error] =
            std::to_chars(buffer.data(), buffer.data() + buffer.size(), value);
        require(error == std::errc{}, capacity);
        append(std::string_view{buffer.data(), end});
    }

    [[nodiscard]] std::size_t size() const noexcept {
        return size_;
    }

private:
    std::size_t limit_;
    std::span<char> output_{};
    std::size_t size_ = 0;
};

void write_value(CsvWriter& csv, const SiEnclosure& enclosure) {
    if (const auto* value = std::get_if<ScaledBinary64>(&enclosure)) {
        require(std::isfinite(value->center) && std::isfinite(value->radius), nonfinite);
        require(value->radius >= 0.0);
        csv.append("binary64,");
        csv.number(value->center);
        csv.append(",,");
        csv.number(value->radius);
        csv.append(",");
        csv.number(value->exponent2);
        return;
    }
    const auto* value = std::get_if<ScaledTwoWord>(&enclosure);
    require(value != nullptr);
    require(std::isfinite(value->high) && std::isfinite(value->low) && std::isfinite(value->radius),
            nonfinite);
    require(value->radius >= 0.0);
    csv.append("two_word,");
    csv.number(value->high);
    csv.append(",");
    csv.number(value->low);
    csv.append(",");
    csv.number(value->radius);
    csv.append(",");
    csv.number(value->exponent2);
}

void write_record(CsvWriter& csv, const ObservationCsvSample& sample,
                  const ObservationCsvColumn& column, const SiEnclosure& value) {
    csv.number(sample.time_seconds.hi);
    csv.append(",");
    csv.number(sample.time_seconds.lo);
    csv.append(",");
    csv.number(sample.checkpoint);
    csv.append(",");
    csv.number(sample.segment);
    csv.append(",");
    csv.number(column.observation.value());
    csv.append(",");
    csv.number(column.component);
    csv.append(",");
    for (const auto label : {column.quantity, column.unit, column.basis, column.reference_point}) {
        csv.label(label);
        csv.append(",");
    }
    if (column.layer)
        csv.number(column.layer->value());
    csv.append(",");
    csv.append(side_name(column.side));
    csv.append(",");
    write_value(csv, value);
    csv.append("\r\n");
}

void write_series(CsvWriter& csv, std::span<const ObservationCsvColumn> columns,
                  std::span<const ObservationCsvSample> samples, const std::stop_token& stop) {
    poll(stop);
    csv.append(header);
    std::optional<std::uint64_t> previous;
    for (const auto& sample : samples) {
        poll(stop);
        require(!previous || *previous < sample.checkpoint);
        require(sample.values.size() == columns.size());
        require(std::isfinite(sample.time_seconds.hi) && std::isfinite(sample.time_seconds.lo),
                nonfinite);
        for (std::size_t i = 0; i < columns.size(); ++i) {
            poll(stop);
            write_record(csv, sample, columns[i], sample.values[i]);
        }
        previous = sample.checkpoint;
    }
    poll(stop);
}

} // namespace

core::Result<core::PlannedBuffer<char>>
encode_observation_csv(std::span<const ObservationCsvColumn> columns,
                       std::span<const ObservationCsvSample> samples, ObservationCsvLimits limits,
                       core::ResourcePlan& resources, const std::stop_token& stop) noexcept {
    return core::guard_entry(operation, invalid, {},
                             [&]() -> core::Result<core::PlannedBuffer<char>> {
                                 try {
                                     poll(stop);
                                     require(columns.size() <= limits.maximum_columns &&
                                                 samples.size() <= limits.maximum_samples,
                                             capacity);
            require(!columns.empty() || samples.empty());
                                     validate_columns(columns, limits, resources, stop);
                                     CsvWriter counter{limits.maximum_output_bytes};
                                     write_series(counter, columns, samples, stop);
                                     auto output = allocate<char>(resources, counter.size());
                                     CsvWriter writer{output.values()};
                                     write_series(writer, columns, samples, stop);
                                     return output;
                                 } catch (const core::Failure& failure) {
                                     return std::unexpected(failure);
                                 }
                             });
}

} // namespace uvlas::io::detail
