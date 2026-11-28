#include "observation_condition.hpp"

#include "input_document_detail.hpp"
#include "source_digest.hpp"

#include <algorithm>
#include <charconv>
#include <limits>

namespace uvlas::io::detail {
namespace {

constexpr auto operation = core::Operation::PrepareModel;

struct Abort {
    core::Failure failure;
};

[[noreturn]] void fail(core::FailureCode code) {
    throw Abort{{code, operation}};
}

template <class T>
T take(core::Result<T> result) {
    if (!result)
        throw Abort{result.error()};
    return std::move(*result);
}

struct Framer {
    Framer(const std::stop_token& cancel, std::uint64_t work, std::uint64_t hash_bytes) noexcept
        : stop(cancel), remaining_work(work), remaining_hash(hash_bytes) {}

    boost::hash2::sha2_256 hash{};
    const std::stop_token& stop;
    std::uint64_t remaining_work;
    std::uint64_t remaining_hash;

    void work(std::uint64_t amount) {
        if (stop.stop_requested())
            fail(core::FailureCode::Cancelled);
        if (amount > remaining_work)
            fail(core::FailureCode::ResourceExceeded);
        remaining_work -= amount;
    }

    void bytes(std::string_view text) {
        work(text.size());
        if (text.size() > remaining_hash)
            fail(core::FailureCode::ResourceExceeded);
        remaining_hash -= text.size();
        if (!text.empty())
            hash.update(text.data(), text.size());
    }

    void frame(std::string_view tag, std::uint64_t count) {
        std::array<char, 24> buffer{};
        const auto end = std::to_chars(buffer.data(), buffer.data() + buffer.size(), count).ptr;
        bytes(tag);
        bytes({buffer.data(), static_cast<std::size_t>(end - buffer.data())});
        bytes(":");
    }

    void text(std::string_view value) {
        frame("s", value.size());
        bytes(value);
    }

    void number(std::string_view token) {
        work(token.size());
        const bool negative = token.front() == '-';
        if (negative)
            token.remove_prefix(1);
        const auto at = token.find_first_of("eE");
        const auto coefficient = token.substr(0, at);
        const auto dot = coefficient.find('.');
        const auto fractional = dot == coefficient.npos ? 0 : coefficient.size() - dot - 1;
        const auto first = coefficient.find_first_not_of("0.");
        if (first == coefficient.npos) {
            frame("d", 1);
            bytes("0");
            return;
        }
        const auto last = coefficient.find_last_not_of("0.");
        const auto trailing = coefficient.size() - last - 1 -
                              static_cast<std::size_t>(dot != coefficient.npos && dot > last);
        std::int64_t exponent = 0;
        if (at != token.npos) {
            auto tail = token.substr(at + 1);
            if (!tail.empty() && tail.front() == '+')
                tail.remove_prefix(1);
            const auto result = std::from_chars(tail.data(), tail.data() + tail.size(), exponent);
            if (result.ec == std::errc::result_out_of_range)
                fail(core::FailureCode::ResourceExceeded);
            if (result.ec != std::errc{} || result.ptr != tail.data() + tail.size())
                fail(core::FailureCode::InvalidInput);
        }
        const auto change =
            static_cast<std::int64_t>(trailing) - static_cast<std::int64_t>(fractional);
        if ((change < 0 && exponent < std::numeric_limits<std::int64_t>::min() - change) ||
            (change > 0 && exponent > std::numeric_limits<std::int64_t>::max() - change))
            fail(core::FailureCode::ResourceExceeded);
        exponent += change;
        std::array<char, 24> encoded{};
        const auto end =
            std::to_chars(encoded.data(), encoded.data() + encoded.size(), exponent).ptr;
        const std::string_view suffix{encoded.data(),
                                      static_cast<std::size_t>(end - encoded.data())};
        const bool split = dot != coefficient.npos && first < dot && dot < last;
        frame("d", last - first + 1 - static_cast<std::size_t>(split) +
                       static_cast<std::size_t>(negative) + 1 + suffix.size());
        if (negative)
            bytes("-");
        if (split) {
            bytes(coefficient.substr(first, dot - first));
            bytes(coefficient.substr(dot + 1, last - dot));
        } else {
            bytes(coefficient.substr(first, last - first + 1));
        }
        bytes("e");
        bytes(suffix);
    }

    template <class Integer>
    void native_number(Integer value) {
        std::array<char, 24> buffer{};
        const auto end = std::to_chars(buffer.data(), buffer.data() + buffer.size(), value).ptr;
        number({buffer.data(), static_cast<std::size_t>(end - buffer.data())});
    }
};

struct PendingValue {
    core::Index node = -1;
    std::optional<std::string_view> key;
};

std::string_view view_text(const InputDocumentData& data, TextRange range) {
    return {data.text.values().data() + range.offset, range.size};
}

bool byte_less(std::string_view a, std::string_view b) {
    return std::lexicographical_compare(
        a.begin(), a.end(), b.begin(), b.end(),
        [](unsigned char lhs, unsigned char rhs) { return lhs < rhs; });
}

} // namespace

core::Result<SourceDigest> observation_condition_digest(const InputNodeView& root,
                                                        ObservationConditionLimits limits,
                                                        core::ResourcePlan& plan,
                                                        const std::stop_token& stop) noexcept {
    return core::guard_entry(
        operation, core::FailureCode::InvalidInput, {}, [&]() -> core::Result<SourceDigest> {
            try {
                Framer frame{stop, limits.maximum_work,
                             std::min(limits.maximum_hash_bytes,
                                      std::numeric_limits<std::uint64_t>::max() / 8)};
                frame.work(1);
                const auto* data = InputDocumentFactory::storage(root);
                if (data == nullptr)
                    fail(core::FailureCode::InvalidInput);
                const auto object = root.object();
                if (!object)
                    fail(core::FailureCode::InvalidInput);
                if (!std::in_range<core::Index>(data->node_count) ||
                    !std::in_range<std::int64_t>(data->text_count))
                    fail(core::FailureCode::ResourceExceeded);
                auto stack = take(plan.allocate<PendingValue>(
                    static_cast<core::Index>(data->node_count), operation));
                std::size_t used = 0;
                constexpr std::array<std::string_view, 8> root_keys{
                    "analysis", "flow",    "initial", "model",
                    "numerics", "restart", "schema",  "units"};
                for (auto key = root_keys.rbegin(); key != root_keys.rend(); ++key) {
                    frame.work(static_cast<std::uint64_t>(object->size()) + 1);
                    const auto child = object->find(*key);
                    if (child) {
                        stack.values()[used++] = {InputDocumentFactory::position(*child), *key};
                    } else if (*key != "initial" && *key != "restart") {
                        fail(core::FailureCode::InvalidInput);
                    }
                }
                frame.bytes("uvlas.observation-condition/1\n");
                frame.frame("o", used);
                while (used != 0) {
                    frame.work(1);
                    const auto item = stack.values()[--used];
                    if (item.key)
                        frame.text(*item.key);
                    const auto& value =
                        data->nodes.values()[static_cast<std::size_t>(item.node)].value;
                    if (std::holds_alternative<std::monostate>(value)) {
                        frame.bytes("n");
                    } else if (const auto* boolean = std::get_if<bool>(&value)) {
                        frame.bytes(*boolean ? "t" : "f");
                    } else if (const auto* text = std::get_if<StringValue>(&value)) {
                        frame.text(view_text(*data, text->text));
                    } else if (const auto* signed_integer = std::get_if<std::int64_t>(&value)) {
                        frame.native_number(*signed_integer);
                    } else if (const auto* unsigned_integer = std::get_if<std::uint64_t>(&value)) {
                        frame.native_number(*unsigned_integer);
                    } else if (const auto* decimal = std::get_if<DecimalValue>(&value)) {
                        frame.number(view_text(*data, decimal->text));
                    } else if (const auto* array = std::get_if<ArrayValue>(&value)) {
                        const auto count = static_cast<std::size_t>(array->children.size);
                        frame.frame("a", count);
                        auto at = used + count;
                        for (auto index = array->children.head; index >= 0;) {
                            frame.work(1);
                            const auto& edge =
                                data->edges.values()[static_cast<std::size_t>(index)];
                            stack.values()[--at] = {edge.child, std::nullopt};
                            index = edge.next;
                        }
                        used += count;
                    } else {
                        const auto& members = std::get<ObjectValue>(value).members;
                        frame.frame("o", static_cast<std::uint64_t>(members.size));
                        const auto first = used;
                        for (auto index = members.head; index >= 0;) {
                            frame.work(1);
                            const auto& member =
                                data->members.values()[static_cast<std::size_t>(index)];
                            stack.values()[used++] = {member.child, view_text(*data, member.key)};
                            index = member.next;
                        }
                        std::ranges::sort(stack.values().subspan(first, used - first),
                                          [&](const PendingValue& a, const PendingValue& b) {
                                              frame.work(1 +
                                                         std::min(a.key->size(), b.key->size()));
                                              return byte_less(*b.key, *a.key);
                                          });
                    }
                }
                frame.work(1);
                return finish_source_digest(frame.hash);
            } catch (const Abort& failure) {
                return std::unexpected(failure.failure);
            }
        });
}

} // namespace uvlas::io::detail
