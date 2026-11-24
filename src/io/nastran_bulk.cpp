#include "nastran_bulk.hpp"

#include <algorithm>
#include <limits>
#include <utility>

namespace uvlas::io::detail {
namespace {

using Code = core::FailureCode;
using Format = NastranBulkFormat;
using Reason = NastranBulkFailure;
using Position = NastranBulkPosition;
using Lexeme = NastranBulkLexeme;
constexpr auto operation = core::Operation::PrepareModel;
constexpr auto no_row = std::numeric_limits<std::size_t>::max();
constexpr auto field_end = std::size_t{72};
constexpr auto line_end = std::size_t{80};

struct Abort {
    core::Failure failure;
    NastranBulkDiagnostic diagnostic;
};

[[noreturn]] void reject(Code code, Reason reason, Position position = {}) {
    throw Abort{{code, operation}, {reason, position}};
}

class Work {
public:
    Work(std::uint64_t limit, const std::stop_token& stop) noexcept : limit_(limit), stop_(stop) {}

    void charge(std::uint64_t amount, Position position = {}) {
        if (stop_.stop_requested())
            reject(Code::Cancelled, Reason::Limit, position);
        if (amount > limit_ - used_)
            reject(Code::ResourceExceeded, Reason::Limit, position);
        used_ += amount;
    }

    [[nodiscard]] std::uint64_t used() const noexcept {
        return used_;
    }

private:
    std::uint64_t limit_;
    const std::stop_token& stop_;
    std::uint64_t used_ = 0;
};

[[nodiscard]] bool letter(char c) noexcept {
    return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z');
}

[[nodiscard]] bool alphanumeric(char c) noexcept {
    return letter(c) || (c >= '0' && c <= '9');
}

[[nodiscard]] bool same_keyword(std::string_view lhs, std::string_view rhs) noexcept {
    return std::ranges::equal(lhs, rhs, [](char actual, char expected) {
        return (actual >= 'a' && actual <= 'z' ? actual - ('a' - 'A') : actual) == expected;
    });
}

[[nodiscard]] std::string_view spelling(std::string_view source, Lexeme token) noexcept {
    return source.substr(token.position.byte, token.size);
}

struct Line {
    std::string_view text;
    std::size_t offset;
    std::size_t number;

    [[nodiscard]] Position at(std::size_t column) const noexcept {
        return {offset + std::min(column, text.size()), number, column + 1};
    }

    [[nodiscard]] Lexeme field(std::size_t first, std::size_t width) const noexcept {
        if (first >= text.size())
            return {at(first), 0};
        const auto raw = text.substr(first, width);
        const auto begin = raw.find_first_not_of(' ');
        if (begin == raw.npos)
            return {at(first), 0};
        const auto last = raw.find_last_not_of(' ');
        return {at(first + begin), last - begin + 1};
    }
};

enum class RowKind : std::uint8_t { Card, Automatic, Named };

struct Row {
    Lexeme physical;
    Lexeme first;
    Lexeme tail;
    Lexeme incoming;
    Lexeme outgoing;
    std::size_t next = no_row;
    std::size_t count = 0;
    Format format = Format::Small;
    RowKind kind = RowKind::Card;
    bool used = false;
};

[[nodiscard]] Lexeme label(const Line& line, std::size_t first) noexcept {
    const auto text = line.field(first, 7);
    const auto value = line.text.substr(text.position.byte - line.offset, text.size);
    return text.position.column == first + 1 && std::ranges::all_of(value, alphanumeric)
               ? text
               : Lexeme{line.at(first), 0};
}

template <class Accept>
void data_fields(const Line& line, Format format, Work& work, Accept accept) {
    const auto produce = [&](Lexeme token) {
        work.charge(1, token.position);
        const auto value = line.text.substr(token.position.byte - line.offset, token.size);
        if (value.contains(' '))
            reject(Code::InvalidInput, Reason::Field, token.position);
        accept(token);
    };
    if (format == Format::Free) {
        for (auto comma = line.text.find(','); comma != line.text.npos;) {
            const auto first = comma + 1;
            comma = line.text.find(',', first);
            const auto end = comma == line.text.npos ? line.text.size() : comma;
            produce(line.field(first, end - first));
        }
        return;
    }
    const auto width = format == Format::Large ? std::size_t{16} : std::size_t{8};
    for (std::size_t first = 8; first < field_end; first += width)
        produce(line.field(first, width));
}

[[nodiscard]] Row describe(const Line& line, Work& work) {
    Row row;
    row.physical = {line.at(0), line.text.size()};
    const auto header = line.text.substr(0, 8);
    const auto separator = header.find_first_of(",=");
    const bool free = separator != header.npos;
    if (free && header[separator] == '=')
        reject(Code::Unsupported, Reason::UnsupportedForm, line.at(separator));
    row.first = line.field(0, free ? separator : 8);
    const auto first = line.text.substr(row.first.position.byte - line.offset, row.first.size);
    if (free) {
        row.format = Format::Free;
        if (first.empty()) {
            if (separator != 0)
                reject(Code::InvalidInput, Reason::Header, row.first.position);
            row.kind = RowKind::Automatic;
        } else if (first.front() == '+' || first.front() == '*' || first.back() == '*') {
            reject(Code::Unsupported, Reason::UnsupportedForm, row.first.position);
        }
    } else {
        row.tail = line.field(field_end, 8);
        row.outgoing = label(line, field_end + 1);
        if (first.empty()) {
            row.kind = RowKind::Automatic;
        } else if (first.front() == '+' || first.front() == '*') {
            row.kind = RowKind::Named;
            row.format = first.front() == '*' ? Format::Large : Format::Small;
            row.incoming = label(line, 1);
            if (row.first.position.column != 1 || row.incoming.size == 0)
                reject(Code::InvalidInput, Reason::Continuation, row.first.position);
        } else if (first.back() == '*') {
            row.format = Format::Large;
            --row.first.size;
        }
    }
    if (row.kind == RowKind::Card) {
        const auto name = first.substr(0, row.first.size);
        if (row.first.position.column != 1 || name.empty() || !letter(name.front()) ||
            !std::ranges::all_of(name, alphanumeric))
            reject(Code::InvalidInput, Reason::Header, row.first.position);
    }
    data_fields(line, row.format, work, [&](Lexeme token) {
        if (row.format == Format::Free) {
            const auto value = line.text.substr(token.position.byte - line.offset, token.size);
            if (value.size() > 1 && (value.front() == '+' || value.front() == '*') &&
                letter(value[1]))
                reject(Code::Unsupported, Reason::UnsupportedForm, token.position);
        }
        ++row.count;
    });
    if (row.format == Format::Free && row.count > 8 &&
        line.text.find_last_not_of(' ') + 1 >= line_end)
        reject(Code::Unsupported, Reason::UnsupportedForm, line.at(line_end - 1));
    return row;
}

template <class Visit>
void scan(std::string_view bytes, NastranBulkLimits limits, Work& work, Visit visit) {
    bool began = false;
    bool ended = false;
    work.charge(bytes.size());
    for (std::size_t offset = 0, number = 1; offset < bytes.size(); ++number) {
        auto newline = bytes.find('\n', offset);
        if (newline == bytes.npos)
            newline = bytes.size();
        const auto raw = bytes.substr(offset, newline - offset);
        if (raw.size() > limits.physical_line_bytes || number > limits.physical_lines)
            reject(Code::ResourceExceeded, Reason::Limit, {offset, number, 1});
        work.charge(raw.size(), {offset, number, 1});
        auto text = raw;
        if (newline < bytes.size() && text.ends_with('\r'))
            text.remove_suffix(1);
        if (const auto null = text.find('\0'); null != text.npos)
            reject(Code::InvalidInput, Reason::SourceByte, {offset + null, number, null + 1});
        text = text.substr(0, text.find('$'));
        if (const auto last = text.find_last_not_of(' '); last != text.npos) {
            const auto body = text.substr(0, last + 1);
            const Line line{text.substr(0, line_end), offset, number};
            if (body.size() > line_end)
                reject(Code::Unsupported, Reason::LineWidth, line.at(line_end));
            for (std::size_t i = 0; i < body.size(); ++i) {
                const auto c = static_cast<unsigned char>(body[i]);
                if (c < 32 || c > 126)
                    reject(Code::Unsupported, Reason::SourceByte, line.at(i));
            }
            if (ended)
                reject(Code::InvalidInput, Reason::Header, line.at(0));

            const auto word_end = body.find_first_of(" ,=");
            const auto word = body.substr(0, word_end);
            if (same_keyword(body, "BEGIN BULK")) {
                if (began)
                    reject(Code::InvalidInput, Reason::Header, line.at(0));
                began = true;
            } else if (same_keyword(body, "ENDDATA")) {
                ended = true;
            } else {
                for (const auto reserved : {"INCLUDE", "SOL", "CEND", "BEGIN", "ENDDATA", "TITLE",
                                            "SUBCASE", "NASTRAN", "ASSIGN"}) {
                    if (same_keyword(word, reserved))
                        reject(Code::Unsupported, Reason::UnsupportedForm, line.at(0));
                }
                began = true;
                visit(describe(line, work));
            }
        }
        offset = newline < bytes.size() ? newline + 1 : newline;
    }
}

struct Counts {
    std::size_t rows = 0;
    std::size_t labels = 0;
    std::size_t cards = 0;
    std::size_t fields = 0;
};

void add(std::size_t& count, std::size_t more, std::uint64_t limit, Position position) {
    if (more > limit - count || more > std::numeric_limits<std::size_t>::max() - count)
        reject(Code::ResourceExceeded, Reason::Limit, position);
    count += more;
}

struct Continuation {
    Lexeme name;
    std::size_t row = 0;
    bool incoming = false;
};

void link(std::span<Row> rows, std::size_t from, std::size_t to) {
    const auto& following = rows[to];
    const auto preceding = rows[from].format;
    if (preceding != following.format &&
        (preceding == Format::Free || following.format == Format::Free))
        reject(Code::Unsupported, Reason::UnsupportedForm, following.first.position);
    if (rows[from].next != no_row)
        reject(Code::InvalidInput, Reason::Continuation, following.first.position);
    rows[from].next = to;
}

void link_named(std::string_view bytes, std::span<Row> rows, std::span<Continuation> labels,
                Work& work) {
    std::ranges::sort(labels, [&](const Continuation& left, const Continuation& right) {
        work.charge(1, left.name.position);
        const auto l = spelling(bytes, left.name);
        const auto r = spelling(bytes, right.name);
        return l < r || (l == r && left.row < right.row);
    });
    while (!labels.empty()) {
        const auto key = spelling(bytes, labels.front().name);
        std::size_t incoming = no_row;
        std::size_t outgoing = no_row;
        std::optional<Position> duplicate_outgoing;
        std::size_t size = 0;
        while (size < labels.size() && spelling(bytes, labels[size].name) == key) {
            const auto& item = labels[size++];
            work.charge(1, item.name.position);
            if (item.incoming) {
                if (incoming != no_row)
                    reject(Code::InvalidInput, Reason::Continuation, item.name.position);
                incoming = item.row;
            } else if (outgoing != no_row) {
                if (!duplicate_outgoing)
                    duplicate_outgoing = item.name.position;
            } else {
                outgoing = item.row;
            }
        }
        if (incoming != no_row) {
            if (duplicate_outgoing)
                reject(Code::InvalidInput, Reason::Continuation, *duplicate_outgoing);
            if (outgoing == no_row)
                reject(Code::InvalidInput, Reason::Continuation, rows[incoming].incoming.position);
            link(rows, outgoing, incoming);
        }
        labels = labels.subspan(size);
    }
}

void connect(std::string_view bytes, std::span<Row> rows, std::span<Continuation> labels,
             Work& work) {
    std::size_t count = 0;
    for (std::size_t i = 0; i < rows.size(); ++i) {
        const auto& row = rows[i];
        work.charge(1, row.first.position);
        if (row.incoming.size != 0)
            labels[count++] = {row.incoming, i, true};
        if (row.outgoing.size != 0)
            labels[count++] = {row.outgoing, i, false};
        if (row.kind == RowKind::Automatic) {
            if (i == 0 || rows[i - 1].tail.size != 0)
                reject(Code::InvalidInput, Reason::Continuation, row.first.position);
            const auto& prior = rows[i - 1];
            // The documented free-field continuation starts after a complete eight-field group.
            // A partial group could require implied padding before the next physical line.
            if (prior.format == Format::Free && row.format == Format::Free && prior.count % 8 != 0)
                reject(Code::Unsupported, Reason::UnsupportedForm, row.first.position);
            link(rows, i - 1, i);
        }
    }
    link_named(bytes, rows, labels, work);
    for (const auto& row : rows) {
        work.charge(1, row.tail.position);
        const auto tail = spelling(bytes, row.tail);
        if (row.next == no_row && !tail.empty() && (tail.front() == '+' || tail.front() == '*'))
            reject(Code::Unsupported, Reason::UnsupportedForm, row.tail.position);
    }
}

void flatten(std::string_view bytes, std::span<Row> rows, std::span<NastranBulkCard> cards,
             std::span<Lexeme> fields, Work& work) {
    std::size_t card_slot = 0;
    std::size_t field_slot = 0;
    for (std::size_t i = 0; i < rows.size(); ++i) {
        const auto& head = rows[i];
        if (head.kind != RowKind::Card)
            continue;
        const auto start = field_slot;
        bool incomplete_large_group = false;
        for (auto next = i; next != no_row; next = rows[next].next) {
            auto& row = rows[next];
            work.charge(1, row.first.position);
            if (row.used)
                reject(Code::InvalidInput, Reason::Continuation, row.first.position);
            row.used = true;
            if (row.format == Format::Large)
                incomplete_large_group = !incomplete_large_group;
            else if (incomplete_large_group)
                reject(Code::Unsupported, Reason::UnsupportedForm, row.first.position);
            const auto pos = row.physical.position;
            data_fields({spelling(bytes, row.physical), pos.byte, pos.line}, row.format, work,
                        [&](Lexeme field) { fields[field_slot++] = field; });
        }
        if (incomplete_large_group)
            reject(Code::InvalidInput, Reason::Continuation, head.first.position);
        cards[card_slot++] = {head.format, head.first, start, field_slot - start};
    }
    for (const auto& row : rows) {
        work.charge(1, row.first.position);
        if (!row.used)
            reject(Code::InvalidInput, Reason::Continuation, row.first.position);
    }
}

template <class T>
void admit(std::size_t count, std::uint64_t& remaining, Work& work) {
    if (!std::in_range<core::Index>(count) || count > remaining / sizeof(T))
        reject(Code::ResourceExceeded, Reason::Limit);
    work.charge(count);
    remaining -= sizeof(T) * count;
}

template <class T>
auto allocate(core::ResourcePlan& plan, std::size_t count) {
    auto result = plan.allocate<T>(static_cast<core::Index>(count), operation);
    if (!result)
        throw Abort{result.error(), {Reason::Limit, {}}};
    return std::move(*result);
}

std::string_view text_of(std::span<const std::byte> bytes) noexcept {
    return {reinterpret_cast<const char*>(bytes.data()), bytes.size()};
}

} // namespace

std::string_view NastranBulkView::text(NastranBulkLexeme lexeme) const noexcept {
    if (lexeme.position.byte > bytes.size() || lexeme.size > bytes.size() - lexeme.position.byte)
        return {};
    return spelling(bytes, lexeme);
}

std::span<const NastranBulkLexeme>
NastranBulkView::data_fields(const NastranBulkCard& card) const noexcept {
    if (card.first_field > fields.size() || card.field_count > fields.size() - card.first_field)
        return {};
    return fields.subspan(card.first_field, card.field_count);
}

NastranBulkDeck::NastranBulkDeck(AcquiredRootSnapshot&& source,
                                 core::PlannedBuffer<NastranBulkCard>&& cards,
                                 core::PlannedBuffer<NastranBulkLexeme>&& fields,
                                 std::uint64_t work) noexcept
    : source_(std::move(source)),
      cards_(std::move(cards)),
      fields_(std::move(fields)),
      work_(work) {}

std::optional<NastranBulkView> NastranBulkDeck::view() const& noexcept {
    const auto view = source_.file.view();
    if (!view)
        return std::nullopt;
    return NastranBulkView{source_.location, view->digest,     text_of(view->bytes),
                           cards_.values(),  fields_.values(), work_};
}

core::Result<NastranBulkDeck>
NastranBulkDeck::prepare(AcquiredRootSnapshot&& snapshot, NastranBulkLimits limits,
                         core::ResourcePlan& plan, const std::stop_token& stop,
                         std::optional<NastranBulkDiagnostic>* diagnostic) noexcept {
    if (diagnostic)
        diagnostic->reset();
    return core::guard_entry(
        operation, Code::InvalidInput, {}, [&]() -> core::Result<NastranBulkDeck> {
            try {
                Work work{limits.work, stop};
                work.charge(0);
                const auto source = snapshot.file.view();
                if (!source || snapshot.location.role != ReferenceRole::NastranDeck)
                    reject(Code::InvalidInput, Reason::Snapshot);
                if (source->bytes.size() > limits.file_bytes)
                    reject(Code::ResourceExceeded, Reason::Limit);
                const auto bytes = text_of(source->bytes);
                Counts counts;
                scan(bytes, limits, work, [&](const Row& row) {
                    add(counts.rows, 1, limits.physical_lines, row.first.position);
                    add(counts.fields, row.count, limits.fields, row.first.position);
                    if (row.kind == RowKind::Card)
                        add(counts.cards, 1, limits.cards, row.first.position);
                    for (const auto label : {row.incoming, row.outgoing}) {
                        if (label.size != 0)
                            add(counts.labels, 1, std::numeric_limits<std::uint64_t>::max(),
                                label.position);
                    }
                });
                if (counts.rows == 0)
                    reject(Code::InvalidInput, Reason::Header);
                auto remaining = plan.remaining_bytes();
                admit<Row>(counts.rows, remaining, work);
                admit<Continuation>(counts.labels, remaining, work);
                admit<NastranBulkCard>(counts.cards, remaining, work);
                admit<Lexeme>(counts.fields, remaining, work);
                auto rows = allocate<Row>(plan, counts.rows);
                auto labels = allocate<Continuation>(plan, counts.labels);
                auto cards = allocate<NastranBulkCard>(plan, counts.cards);
                auto fields = allocate<Lexeme>(plan, counts.fields);
                std::size_t i = 0;
                scan(bytes, limits, work, [&](const Row& row) { rows.values()[i++] = row; });
                connect(bytes, rows.values(), labels.values(), work);
                flatten(bytes, rows.values(), cards.values(), fields.values(), work);
                work.charge(0);
                return NastranBulkDeck{std::move(snapshot), std::move(cards), std::move(fields),
                                       work.used()};
            } catch (const Abort& error) {
                if (diagnostic && error.failure.code != Code::Cancelled)
                    *diagnostic = error.diagnostic;
                return std::unexpected(error.failure);
            }
        });
}

} // namespace uvlas::io::detail
