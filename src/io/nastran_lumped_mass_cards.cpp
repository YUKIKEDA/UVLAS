#include "nastran_lumped_mass_cards.hpp"

#include <algorithm>
#include <limits>
#include <new>
#include <utility>

namespace uvlas::io::detail {
namespace {

using Code = core::FailureCode;
using Reason = NastranLumpedMassFailure;
using Lexeme = NastranBulkLexeme;
using Position = NastranBulkPosition;
constexpr auto operation = core::Operation::PrepareModel;

struct Abort {
    core::Failure failure;
    NastranLumpedMassDiagnostic diagnostic;
};

[[noreturn]] void reject(Reason reason, Position position = {}, Code code = Code::InvalidInput) {
    throw Abort{{code, operation}, {reason, position}};
}

struct Work {
    std::uint64_t limit;
    const std::stop_token& stop;
    std::uint64_t used = 0;

    void visit(Position position = {}) {
        if (stop.stop_requested())
            reject(Reason::Limit, position, Code::Cancelled);
        if (used == limit)
            reject(Reason::Limit, position, Code::ResourceExceeded);
        ++used;
    }
};

[[nodiscard]] bool is_conm2(std::string_view text) noexcept {
    return std::ranges::equal(text, std::string_view{"CONM2"}, [](char actual, char expected) {
        return (actual >= 'a' && actual <= 'z' ? actual - ('a' - 'A') : actual) == expected;
    });
}

template <class T>
[[nodiscard]] auto allocate(core::ResourcePlan& memory, std::size_t size) {
    if (!std::in_range<core::Index>(size))
        reject(Reason::Limit, {}, Code::ResourceExceeded);
    auto result = memory.allocate<T>(static_cast<core::Index>(size), operation);
    if (!result)
        throw Abort{result.error(), {Reason::Limit, {}}};
    return std::move(*result);
}

struct IndexedId {
    std::uint64_t id = 0;
    std::size_t ordinal = 0;
};

template <class T>
[[nodiscard]] auto make_index(std::span<const T> rows, Work& work, core::ResourcePlan& memory) {
    auto result = allocate<IndexedId>(memory, rows.size());
    const auto index = result.values();
    for (std::size_t i = 0; i < rows.size(); ++i) {
        work.visit(rows[i].id_field.position);
        index[i] = {rows[i].id.value(), i};
    }
    std::ranges::sort(
        index,
        [&](std::uint64_t lhs, std::uint64_t rhs) {
            work.visit();
            return lhs < rhs;
        },
        &IndexedId::id);
    return result;
}

[[nodiscard]] std::size_t find_id(std::span<const IndexedId> index, std::uint64_t id, Lexeme field,
                                  Work& work, Reason missing) {
    const auto found = std::ranges::lower_bound(
        index, id,
        [&](std::uint64_t lhs, std::uint64_t rhs) {
            work.visit(field.position);
            return lhs < rhs;
        },
        &IndexedId::id);
    if (found == index.end() || found->id != id)
        reject(missing, field.position);
    return found->ordinal;
}

struct Fields {
    const NastranBulkView& bulk;
    const NastranBulkCard& card;
    std::span<const Lexeme> values;

    [[nodiscard]] std::optional<Lexeme> at(std::size_t index) const noexcept {
        return index < values.size() ? std::optional{values[index]} : std::nullopt;
    }

    [[nodiscard]] Lexeme required(std::size_t index) const {
        const auto field = at(index);
        if (!field || field->size == 0)
            reject(Reason::MissingField, field ? field->position : card.name.position);
        return *field;
    }

    [[nodiscard]] std::int64_t integer(Lexeme field) const {
        if (card.format == NastranBulkFormat::Free && field.size > 8)
            reject(Reason::InvalidField, field.position);
        const auto value = bulk.integer(field);
        if (!value)
            reject(Reason::InvalidField, field.position);
        return *value;
    }

    [[nodiscard]] std::uint64_t positive(Lexeme field) const {
        const auto value = integer(field);
        if (value <= 0)
            reject(Reason::InvalidField, field.position);
        return static_cast<std::uint64_t>(value);
    }

    [[nodiscard]] NastranMassReal real(std::size_t index, Work& work) const {
        const auto field = at(index);
        work.visit(field ? field->position : card.name.position);
        if (!field || field->size == 0)
            return {{}, field};
        const auto value = bulk.real(*field);
        if (!value)
            reject(Reason::InvalidField, field->position);
        return {*value, field};
    }

    void verify_reserved(Work& work) const {
        const auto unused = at(7);
        if (unused && unused->size != 0)
            reject(Reason::ExtraField, unused->position);
        if (values.size() <= 14)
            return;
        for (const auto field : values.subspan(14)) {
            work.visit(field.position);
            if (field.size != 0)
                reject(Reason::ExtraField, field.position);
        }
    }
};

[[nodiscard]] NastranMassFrame mass_frame(const Fields& input, std::span<const IndexedId> frames,
                                          Work& work) {
    const auto field = input.at(2);
    work.visit(field ? field->position : input.card.name.position);
    if (!field || field->size == 0)
        return NastranBasicMassOffset{};
    const auto cid = input.integer(*field);
    if (cid < -1)
        reject(Reason::InvalidField, field->position);
    if (cid == -1)
        return NastranBasicMassCentre{};
    if (cid == 0)
        return NastranBasicMassOffset{};
    const auto id = static_cast<std::uint64_t>(cid);
    return NastranFramedMassOffset{NastranCoordinateId{id},
                                   find_id(frames, id, *field, work, Reason::UnknownFrame)};
}

[[nodiscard]] NastranLumpedMass read_one(const Fields& input, std::span<const IndexedId> grids,
                                         std::span<const IndexedId> frames, Work& work) {
    const auto id = input.required(0);
    const auto grid = input.required(1);
    static_cast<void>(input.required(3)); // M is required; positions and inertia may be blank.
    input.verify_reserved(work);
    const auto gid = input.positive(grid);
    const auto grid_ordinal = find_id(grids, gid, grid, work, Reason::UnknownGrid);
    return {NastranElementId{input.positive(id)},
            NastranGridId{gid},
            grid_ordinal,
            mass_frame(input, frames, work),
            input.card.name,
            id,
            grid,
            input.at(2),
            input.real(3, work),
            {input.real(4, work), input.real(5, work), input.real(6, work)},
            {input.real(8, work), input.real(9, work), input.real(10, work), input.real(11, work),
             input.real(12, work), input.real(13, work)}};
}

} // namespace

NastranLumpedMassCards::NastranLumpedMassCards(NastranLumpedMassSource source,
                                               core::PlannedBuffer<NastranLumpedMass>&& masses,
                                               std::uint64_t work) noexcept
    : source_(source), masses_(std::move(masses)), work_(work) {}

NastranLumpedMassCards::NastranLumpedMassCards(NastranLumpedMassCards&& other) noexcept
    : source_(std::exchange(other.source_, std::nullopt)),
      masses_(std::move(other.masses_)),
      work_(other.work_) {}

core::Result<NastranLumpedMassCards>
NastranLumpedMassCards::prepare(const NastranShellCards& shells, NastranLumpedMassLimits limits,
                                core::ResourcePlan& memory, const std::stop_token& stop,
                                std::optional<NastranLumpedMassDiagnostic>* diagnostic) noexcept {
    if (diagnostic != nullptr)
        diagnostic->reset();
    try {
        Work work{limits.work, stop};
        work.visit();
        const auto held = shells.view();
        if (!held)
            reject(Reason::Input);
        const auto& input = held->input;
        const auto& bulk = input.bulk;
        std::size_t count = 0;
        for (const auto& card : bulk.cards) {
            work.visit(card.name.position);
            if (!is_conm2(bulk.text(card.name)))
                continue;
            if (count >= limits.records || count == std::numeric_limits<std::size_t>::max())
                reject(Reason::Limit, card.name.position, Code::ResourceExceeded);
            ++count;
        }
        auto result = allocate<NastranLumpedMass>(memory, count);
        if (count != 0) {
            const auto grids = make_index(input.grids, work, memory);
            const auto frames = make_index(input.coordinates, work, memory);
            std::size_t slot = 0;
            for (const auto& card : bulk.cards) {
                work.visit(card.name.position);
                if (is_conm2(bulk.text(card.name)))
                    result.values()[slot++] = read_one({bulk, card, bulk.data_fields(card)},
                                                       grids.values(), frames.values(), work);
            }
        }
        return NastranLumpedMassCards{{bulk.location, bulk.digest}, std::move(result), work.used};
    } catch (const Abort& error) {
        if (diagnostic != nullptr && error.failure.code != Code::Cancelled)
            *diagnostic = error.diagnostic;
        return std::unexpected(error.failure);
    } catch (const std::bad_alloc&) {
        return std::unexpected(core::Failure{Code::ResourceExceeded, operation});
    } catch (...) {
        return std::unexpected(core::Failure{Code::InvalidInput, operation});
    }
}

std::optional<NastranLumpedMassView> NastranLumpedMassCards::view() const& noexcept {
    if (!source_)
        return std::nullopt;
    return NastranLumpedMassView{*source_, masses_.values(), work_};
}

} // namespace uvlas::io::detail
