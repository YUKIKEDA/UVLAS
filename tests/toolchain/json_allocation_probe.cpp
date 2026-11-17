// N4-03/N4-08: isolate allocation failures in the fixed JSON SAX dependency.
// This executable never links its process-local allocator into a product target.
#include <nlohmann/json.hpp>

#include <array>
#include <bit>
#include <charconv>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <limits>
#include <new>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

namespace {
constexpr int caught_allocation = 42;
constexpr int unexpected_exception = 43;
constexpr int terminated = 77;
std::int64_t remaining = -1;
std::uint64_t allocation_attempts = 0;

void* allocate(std::size_t size) {
    if (remaining != -1) {
        ++allocation_attempts;
        if (remaining == 0) {
            // Reporting/unwinding must not fail recursively.
            remaining = -1;
            throw std::bad_alloc{};
        }
        if (remaining > 0)
            --remaining;
    }
    // NOLINTNEXTLINE(misc-const-correctness): operator new must return writable pointee storage.
    void* const value = std::malloc(size == 0 ? 1 : size);
    if (value == nullptr)
        throw std::bad_alloc{};
    return value;
}

using Json = nlohmann::json;
struct Events {
    std::array<char, 4096> bytes{};
    std::size_t used = 0;
    std::uint64_t count = 0;
    bool complete = true;

    void append(std::string_view text) noexcept {
        if (text.size() > bytes.size() - used) {
            complete = false;
            return;
        }
        for (const char value : text)
            bytes[used++] = value;
    }
    template <class Integer>
    void integer(Integer value) noexcept {
        std::array<char, 32> text{};
        const auto result = std::to_chars(text.data(), text.data() + text.size(), value);
        if (result.ec != std::errc{}) {
            complete = false;
            return;
        }
        append({text.data(), static_cast<std::size_t>(result.ptr - text.data())});
    }
    void text(char kind, std::string_view value) noexcept {
        append({&kind, 1});
        integer(value.size());
        append(":");
        append(value);
        append(";");
        ++count;
    }
    void marker(std::string_view value) noexcept {
        append(value);
        ++count;
    }
};

// Record actual event kinds, integer values, float bits/raw token, decoded UTF-8
// bytes and ordering without allocating in a user SAX callback.
struct Sax final : nlohmann::json_sax<Json> {
    Events events;
    int error_id = 0;
    std::size_t error_byte = 0;
    bool error_prefix = false;

    bool null() override {
        events.marker("N;");
        return events.complete;
    }
    bool boolean(bool value) override {
        events.marker(value ? "T;" : "F;");
        return events.complete;
    }
    bool number_integer(Json::number_integer_t value) override {
        events.append("I");
        events.integer(value);
        events.marker(";");
        return events.complete;
    }
    bool number_unsigned(Json::number_unsigned_t value) override {
        events.append("U");
        events.integer(value);
        events.marker(";");
        return events.complete;
    }
    bool number_float(Json::number_float_t value, const Json::string_t& raw) override {
        constexpr std::string_view hex = "0123456789abcdef";
        const auto bits = std::bit_cast<std::uint64_t>(value);
        events.append("D");
        for (unsigned nibble = 16; nibble != 0; --nibble) {
            const char digit = hex[static_cast<std::size_t>((bits >> (4U * (nibble - 1))) & 15U)];
            events.append({&digit, 1});
        }
        events.append(":");
        events.integer(raw.size());
        events.append(":");
        events.append(raw);
        events.marker(";");
        return events.complete;
    }
    bool string(Json::string_t& value) override {
        events.text('S', value);
        return events.complete;
    }
    bool binary(Json::binary_t&) override {
        // A JSON-text parser must never manufacture a binary event.
        events.complete = false;
        return false;
    }
    bool start_object(std::size_t) override {
        events.marker("O;");
        return events.complete;
    }
    bool key(Json::string_t& value) override {
        events.text('K', value);
        return events.complete;
    }
    bool end_object() override {
        events.marker("C;");
        return events.complete;
    }
    bool start_array(std::size_t) override {
        events.marker("A;");
        return events.complete;
    }
    bool end_array() override {
        events.marker("Z;");
        return events.complete;
    }
    bool parse_error(std::size_t position, const std::string&,
                     const Json::exception& error) override {
        error_id = error.id;
        error_byte = position;
        error_prefix = std::strncmp(error.what(), "[json.exception.", 16) == 0;
        return false;
    }
};

struct Fixture {
    std::array<char, 1024> bytes{};
    std::size_t size = 0;
    bool valid = true;
    void append(std::string_view value) noexcept {
        if (value.size() > bytes.size() - size) {
            valid = false;
            return;
        }
        for (const char byte : value)
            bytes[size++] = byte;
    }
};
Fixture fixture(std::string_view name) noexcept {
    Fixture result;
    if (name == "all-kinds") {
        result.append(
            R"({"a":[null,true,false,-1,0,9223372036854775807,18446744073709551615,1.5,"a\u0000b"],"b":{"k":"日本語"}})");
    } else if (name == "surrogate") {
        result.append(R"(["\uD800"])");
    } else if (name == "range") {
        result.append("[1e400]");
    } else if (name == "deep") {
        for (int level = 0; level < 80; ++level)
            result.append("[");
        result.append("0");
        for (int level = 0; level < 80; ++level)
            result.append("]");
    } else if (name == "long-token") {
        result.append("[\"");
        for (int byte = 0; byte < 600; ++byte)
            result.append("x");
        result.append("\"] trailing");
    } else {
        result.valid = false;
    }
    return result;
}

constexpr bool debug_proxies =
#if defined(_MSC_VER) && _ITERATOR_DEBUG_LEVEL > 0
    true;
#else
    false;
#endif

void report(bool accepted, const Sax& sax) noexcept {
    constexpr std::string_view hex = "0123456789abcdef";
    std::printf("{\"schema\":1,\"new_count\":%llu,\"accepted\":%s,\"events\":%llu,"
                "\"error_id\":%d,\"error_byte\":%llu,\"error_prefix\":%s,"
                "\"complete\":%s,\"debug_proxies\":%s,\"events_hex\":\"",
                static_cast<unsigned long long>(allocation_attempts), accepted ? "true" : "false",
                static_cast<unsigned long long>(sax.events.count), sax.error_id,
                static_cast<unsigned long long>(sax.error_byte),
                sax.error_prefix ? "true" : "false", sax.events.complete ? "true" : "false",
                debug_proxies ? "true" : "false");
    for (std::size_t i = 0; i < sax.events.used; ++i) {
        const auto byte = static_cast<unsigned char>(sax.events.bytes[i]);
        std::putchar(hex[byte >> 4U]);
        std::putchar(hex[byte & 15U]);
    }
    std::puts("\"}");
}
} // namespace

void* operator new(std::size_t size) {
    return allocate(size);
}
void* operator new[](std::size_t size) {
    return allocate(size);
}
void operator delete(void* value) noexcept {
    std::free(value);
}
void operator delete[](void* value) noexcept {
    std::free(value);
}
void operator delete(void* value, std::size_t) noexcept {
    std::free(value);
}
void operator delete[](void* value, std::size_t) noexcept {
    std::free(value);
}

int main(int argc, char** argv) {
    if (argc != 2 && argc != 3)
        return 2;
    const std::string_view name{argv[1]};
    std::int64_t fault = -2; // Count allocations without injecting a failure.
    if (argc == 3) {
        const std::string_view text{argv[2]};
        const auto parsed = std::from_chars(text.data(), text.data() + text.size(), fault);
        if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size() || fault < 0)
            return 2;
    }
    const bool negative = name == "default-vector" || name == "explicit-bool" ||
                          name == "explicit-char" || name == "explicit-string";
    const auto input = fixture(name);
    if (!negative && !input.valid)
        return 2;
    std::set_terminate([] { std::_Exit(terminated); });
    Sax sax;
    allocation_attempts = 0;
    remaining = fault;
    try {
        bool accepted = true;
        if (name == "default-vector") {
            const std::vector<bool> value;
            accepted = value.empty();
        } else if (name == "explicit-bool") {
            const std::vector<bool> value(std::size_t{0});
            accepted = value.empty();
        } else if (name == "explicit-char") {
            const std::vector<char> value(std::size_t{0});
            accepted = value.empty();
        } else if (name == "explicit-string") {
            const std::string value(std::size_t{0}, char{});
            accepted = value.empty();
        } else {
            // Instantiate JSON text only: binary format readers are outside this
            // patch's call graph. The strict fixed parser still owns all grammar.
            auto adapter = nlohmann::detail::input_adapter(input.bytes.data(),
                                                           input.bytes.data() + input.size);
            nlohmann::detail::parser<Json, decltype(adapter)> parser(decltype(adapter){adapter},
                                                                     nullptr, true, false);
            accepted = parser.sax_parse(&sax, true);
        }
        remaining = -1;
        report(accepted, sax);
        return 0;
    } catch (const std::bad_alloc&) {
        remaining = -1;
        return caught_allocation;
    } catch (...) {
        remaining = -1;
        return unexpected_exception;
    }
}
