#include "resolved_case_document.hpp"

#include "si_case_detail.hpp"
#include "source_digest.hpp"

#include <algorithm>
#include <bit>
#include <charconv>
#include <limits>
#include <type_traits>
#include <utility>

namespace uvlas::io::detail {
namespace {

using Code = core::FailureCode;
constexpr auto operation = core::Operation::WriteResults;
constexpr std::string_view hex_digits = "0123456789abcdef";

[[noreturn]] void fail(Code code) {
    throw core::Failure{code, operation};
}

struct DecimalText {
    template <class Integer>
    explicit DecimalText(Integer value) {
        const auto [end, ec] = std::to_chars(text.data(), text.data() + text.size(), value);
        if (ec != std::errc{}) {
            fail(Code::InvalidInput);
        }
        size = static_cast<std::size_t>(end - text.data());
    }

    [[nodiscard]] std::string_view view() const noexcept {
        return {text.data(), size};
    }
    std::array<char, 24> text{};
    std::size_t size = 0;
};

[[nodiscard]] std::string_view json_escape(char c, std::array<char, 6>& scratch) noexcept {
    if (c == '\\' || c == '"') {
        scratch[0] = '\\';
        scratch[1] = c;
        return {scratch.data(), 2};
    }
    const auto byte = static_cast<unsigned char>(c);
    if (byte < 0x20U) {
        scratch = {'\\', 'u', '0', '0', hex_digits[byte >> 4], hex_digits[byte & 15U]};
        return {scratch.data(), scratch.size()};
    }
    scratch[0] = c;
    return {scratch.data(), 1};
}

class Output {
public:
    Output(std::uint64_t bytes, std::uint64_t work, const std::stop_token& stop,
           std::span<char> storage = {}) noexcept
        : maximum_bytes_(bytes), maximum_work_(work), stop_(stop), storage_(storage) {}

    [[nodiscard]] std::uint64_t size() const noexcept {
        return size_;
    }
    [[nodiscard]] std::uint64_t work() const noexcept {
        return work_;
    }
    [[nodiscard]] std::uint64_t remaining_bytes() const noexcept {
        return maximum_bytes_ - size_;
    }

    void step(std::uint64_t work) {
        if (stop_.stop_requested()) {
            fail(Code::Cancelled);
        }
        if (work > maximum_work_ - work_) {
            fail(Code::ResourceExceeded);
        }
        work_ += work;
    }

    [[nodiscard]] std::span<char> reserve(std::uint64_t count) {
        if (count > remaining_bytes()) {
            fail(Code::ResourceExceeded);
        }
        const auto start = size_;
        step(count);
        size_ += count;
        if (storage_.empty()) {
            return {};
        }
        return storage_.subspan(static_cast<std::size_t>(start), static_cast<std::size_t>(count));
    }

    void append(std::string_view text) {
        if (const auto bytes = reserve(text.size()); !bytes.empty()) {
            std::ranges::copy(text, bytes.begin());
        }
    }

    void quoted(std::string_view text) {
        append("\"");
        std::array<char, 6> scratch{};
        for (const char c : text) {
            step(1);
            append(json_escape(c, scratch));
        }
        append("\"");
    }

    template <class Integer>
    void number(Integer value) {
        const DecimalText text{value};
        append(text.view());
    }

    template <class Integer>
    void identifier(Integer value) {
        append("\"");
        number(value);
        append("\"");
    }

    void boolean(bool value) {
        append(value ? "true" : "false");
    }

    template <std::size_t N>
    void hex(const std::array<std::byte, N>& value) {
        std::array<char, 2 * N> text{};
        for (std::size_t i = 0; i < N; ++i) {
            const auto byte = std::to_integer<unsigned>(value[i]);
            text[2 * i] = hex_digits[byte >> 4];
            text[2 * i + 1] = hex_digits[byte & 15U];
        }
        append("\"");
        append({text.data(), text.size()});
        append("\"");
    }

    void word(double value) {
        const auto word = std::bit_cast<std::uint64_t>(value);
        std::array<std::byte, 8> bytes{};
        for (std::size_t i = 0; i < bytes.size(); ++i) {
            bytes[i] = static_cast<std::byte>((word >> (56U - 8U * i)) & 0xffU);
        }
        hex(bytes);
    }

private:
    std::uint64_t maximum_bytes_, maximum_work_;
    const std::stop_token& stop_;
    std::span<char> storage_;
    std::uint64_t size_ = 0, work_ = 0;
};

template <class Emit>
void reversed_pointer(Output& out, const SchemaData& schema, const InputDocumentData& document,
                      std::uint64_t node, Emit emit) {
    const auto nodes = schema.nodes.values();
    std::array<char, 6> escape{};
    while (const auto parent = nodes[static_cast<std::size_t>(node)].parent) {
        out.step(1);
        const auto& current = nodes[static_cast<std::size_t>(node)];
        if (const auto member = current.parent_member) {
            const auto& key = document.members.values()[static_cast<std::size_t>(*member)].key;
            const auto source = document.text.values().subspan(key.offset, key.size);
            for (auto position = source.rbegin(); position != source.rend(); ++position) {
                out.step(1);
                if (*position == '~') {
                    emit("~0");
                } else if (*position == '/') {
                    emit("~1");
                } else {
                    emit(json_escape(*position, escape));
                }
            }
        } else if (const auto index = current.parent_array_index) {
            const DecimalText text{*index};
            emit(text.view());
        } else {
            fail(Code::InvalidInput);
        }
        emit("/");
        node = *parent;
    }
}

void pointer(Output& out, const SchemaData& schema, const InputDocumentData& document,
             std::uint64_t node) {
    out.append("\"");
    std::uint64_t count = 0;
    reversed_pointer(out, schema, document, node, [&](std::string_view part) {
        if (part.size() > out.remaining_bytes() - count) {
            fail(Code::ResourceExceeded);
        }
        count += part.size();
    });
    const auto destination = out.reserve(count);
    reversed_pointer(out, schema, document, node, [&](std::string_view part) {
        count -= part.size();
        if (!destination.empty()) {
            std::ranges::copy(part, destination.subspan(static_cast<std::size_t>(count)).begin());
        }
    });
    out.append("\"");
}

template <class Unit, std::size_t N>
void unit(Output& out, Unit value, const std::array<std::string_view, N>& names) {
    const auto key = static_cast<std::size_t>(value);
    if (key >= names.size()) {
        fail(Code::InvalidInput);
    }
    out.quoted(names[key]);
}

void units(Output& out, SiUnits value) {
    out.append("{\"length\":");
    unit(out, value.length, std::to_array<std::string_view>({"m", "mm", "cm", "km", "in", "ft"}));
    out.append(",\"mass\":");
    unit(out, value.mass, std::to_array<std::string_view>({"kg", "g", "t", "lbm"}));
    out.append(",\"time\":");
    unit(out, value.time, std::to_array<std::string_view>({"s", "ms", "min"}));
    out.append(",\"angle\":");
    unit(out, value.angle, std::to_array<std::string_view>({"rad", "deg"}));
    out.append("}");
}

void dimensions(Output& out, SchemaDimension value) {
    out.append("[");
    out.number(value.exponents[0]);
    for (const auto exponent : std::span{value.exponents}.subspan(1)) {
        out.append(",");
        out.number(exponent);
    }
    out.append("]");
}

void decimal(Output& out, SiDecimalView value) {
    out.append("{\"sign\":");
    out.number(value.sign);
    out.append(",\"digits\":");
    out.quoted(value.digits);
    out.append(",\"exponent10\":");
    out.identifier(value.exponent10);
    out.append("}");
}

void original(Output& out, const SiRecipeView& recipe) {
    if (const auto* value = std::get_if<SiDecimalView>(&recipe.value)) {
        out.append("{\"decimal\":");
        decimal(out, *value);
        out.append(",\"token\":");
        out.quoted(recipe.original);
    } else {
        const auto& rational = std::get<SiRationalView>(recipe.value);
        out.append("{\"numerator\":");
        decimal(out, rational.numerator);
        out.append(",\"denominator\":");
        decimal(out, rational.denominator);
    }
    out.append("}");
}

void requirement_value(Output& out, const SiExactInput& input) {
    if (const auto* token = std::get_if<SiDecimalInput>(&input)) {
        out.quoted(token->token);
    } else {
        const auto& rational = std::get<SiRationalInput>(input);
        out.append("{\"numerator\":");
        out.quoted(rational.numerator);
        out.append(",\"denominator\":");
        out.quoted(rational.denominator);
        out.append("}");
    }
}

void enclosure(Output& out, const SiEnclosure& value) {
    std::visit(
        [&](const auto& words) {
            if constexpr (std::is_same_v<std::remove_cvref_t<decltype(words)>, ScaledBinary64>) {
                out.append("{\"representation\":\"Binary64\",\"center_bits\":");
                out.word(words.center);
            } else {
                out.append("{\"representation\":\"TwoWord\",\"high_bits\":");
                out.word(words.high);
                out.append(",\"low_bits\":");
                out.word(words.low);
            }
            out.append(",\"radius_bits\":");
            out.word(words.radius);
            out.append(",\"exponent2\":");
            out.identifier(words.exponent2);
            out.append("}");
        },
        value);
}

void write_row(Output& out, const SiData& data, const SchemaData& schema,
               const InputDocumentData& document, std::size_t row) {
    const auto& location = data.entries.values()[row];
    const auto recipe = si_recipe(document, schema.nodes.values(), schema.digits.values(),
                                  location.node, data.units);
    const auto request =
        si_requirement_view(data.requirements.values()[row], data.request_text.values());
    const auto& record = schema_records()[location.record];
    out.append("{\"path\":");
    pointer(out, schema, document, location.node);
    out.append(",\"record\":");
    out.quoted(record.name);
    out.append(",\"field\":");
    out.quoted(record.fields[location.field].name);
    out.append(",\"ledger_origin\":");
    out.identifier(row);
    out.append(",\"original\":");
    original(out, recipe);
    out.append(",\"dimension\":");
    dimensions(out, recipe.dimensions);
    out.append(",\"pi_power\":");
    out.number(recipe.pi_power);
    out.append(",\"enclosure\":");
    enclosure(out, data.converted.values()[row]);
    out.append(",\"requirements\":{\"si_scale\":");
    requirement_value(out, request.scale.value);
    out.append(",\"absolute\":");
    requirement_value(out, request.absolute);
    out.append(",\"relative\":");
    requirement_value(out, request.relative);
    out.append("}}");
}

void write_document(Output& out, const SiData& data, const SchemaData& schema,
                    const SchemaIdentity& identity) {
    const auto& document = schema.snapshot.values()[0];
    out.append("{\"source\":{\"id\":");
    out.identifier(identity.source.source.value());
    out.append(",\"byte_count\":");
    out.identifier(identity.source.byte_count);
    out.append(",\"sha256\":");
    out.hex(identity.source.digest.bytes);
    out.append(",\"contract_revision\":");
    out.number(identity.contract_revision);
    out.append(",\"representation_revision\":");
    out.number(identity.representation_revision);
    out.append("},\"original_case_utf8\":");
    const auto original_bytes = document.source.values();
    out.quoted({original_bytes.data(), original_bytes.size()});
    out.append(",\"original_units\":");
    units(out, data.units);
    out.append(",\"structural_defaults\":{\"provenance_missing\":");
    out.boolean(data.defaults.provenance_missing);
    out.append(",\"initial_missing\":");
    out.boolean(data.defaults.initial_missing);
    out.append(",\"restart_missing\":");
    out.boolean(data.defaults.restart_missing);
    if (data.defaults.provenance_missing) {
        out.append(
            ",\"applied_provenance\":{\"description\":\"\",\"references\":[],\"assumptions\":[]}");
    }
    out.append("},\"quantities\":[");
    for (std::size_t row = 0; row < data.entries.size(); ++row) {
        if (row != 0) {
            out.append(",");
        }
        write_row(out, data, schema, document, row);
    }
    out.append("]}\n");
}

[[nodiscard]] SourceDigest hash(std::span<const char> text, const std::stop_token& stop) {
    constexpr auto batch_bytes = std::size_t{64} * 1024;
    boost::hash2::sha2_256 digest;
    while (!text.empty()) {
        if (stop.stop_requested()) {
            fail(Code::Cancelled);
        }
        const auto count = std::min(batch_bytes, text.size());
        digest.update(text.data(), count);
        text = text.subspan(count);
    }
    if (stop.stop_requested()) {
        fail(Code::Cancelled);
    }
    return finish_source_digest(digest);
}

} // namespace

ResolvedCaseDocument::ResolvedCaseDocument(core::PlannedBuffer<char> text, SourceIdentity source,
                                           SourceDigest digest) noexcept
    : text_(std::move(text)), source_(source), digest_(digest) {}

std::string_view ResolvedCaseDocument::text() const& noexcept {
    const auto text = text_.values();
    return text.empty() ? std::string_view{} : std::string_view{text.data(), text.size()};
}

std::optional<ResultProvenance> ResolvedCaseDocument::provenance() const noexcept {
    if (text_.size() == 0) {
        return std::nullopt;
    }
    return ResultProvenance{source_, digest_};
}

core::Result<ResolvedCaseDocument>
write_resolved_case_document(const SiCase& input, ResolvedCaseDocumentLimits limits,
                             core::ResourcePlan& plan, const std::stop_token& stop) noexcept {
    return core::guard_entry(
        operation, Code::InvalidInput, {}, [&]() -> core::Result<ResolvedCaseDocument> {
            try {
                if (stop.stop_requested()) {
                    fail(Code::Cancelled);
                }
                const auto* data = SiFactory::storage(input);
                if (!data || !data->source) {
                    fail(Code::InvalidInput);
                }
                const auto* schema = SiFactory::storage(*data->source);
                const auto identity = data->source->identity();
                if (!identity || !schema || schema->snapshot.size() != 1) {
                    fail(Code::InvalidInput);
                }
                if (data->entries.size() > limits.maximum_quantities) {
                    fail(Code::ResourceExceeded);
                }
                const auto bytes =
                    std::min({limits.maximum_output_bytes, plan.remaining_bytes(),
                              static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max()),
                              static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())});
                Output measure{bytes, limits.maximum_work, stop};
                write_document(measure, *data, *schema, *identity);
                const auto size = measure.size();
                if (size > limits.maximum_work ||
                    measure.work() > (limits.maximum_work - size) / 2U) {
                    fail(Code::ResourceExceeded);
                }
                if (stop.stop_requested()) {
                    fail(Code::Cancelled);
                }
                auto text = plan.allocate<char>(static_cast<core::Index>(size), operation);
                if (!text) {
                    return std::unexpected(text.error());
                }
                Output output{size, measure.work(), stop, text->values()};
                write_document(output, *data, *schema, *identity);
                if (size != output.size()) {
                    fail(Code::InvalidInput);
                }
                const auto digest = hash(text->values(), stop);
                return ResolvedCaseDocument{std::move(*text), identity->source, digest};
            } catch (const core::Failure& error) {
                return std::unexpected(error);
            }
        });
}

core::Result<ResolvedResultDocuments>
write_resolved_result_documents(const SiCase& input, const ResultStatus& status,
                                ResolvedCaseDocumentLimits limits,
                                std::uint64_t maximum_result_bytes, core::ResourcePlan& plan,
                                const std::stop_token& stop) noexcept {
    auto resolved = write_resolved_case_document(input, limits, plan, stop);
    if (!resolved) {
        return std::unexpected(resolved.error());
    }
    const auto provenance = resolved->provenance();
    if (!provenance) {
        return std::unexpected(core::Failure{Code::InvalidInput, operation});
    }
    auto result = write_result_documents(status, *provenance, maximum_result_bytes, plan, stop);
    if (!result) {
        return std::unexpected(result.error());
    }
    return ResolvedResultDocuments{std::move(*resolved), std::move(*result)};
}

} // namespace uvlas::io::detail
