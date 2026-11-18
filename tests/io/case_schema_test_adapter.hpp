#pragma once

#include "case_schema_corpus.hpp"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace uvlas::test::schema_adapter {

// This fixture editor treats every scalar as its original wire token. In particular,
// no JSON number is parsed as double or reserialized through a binary value.
struct JsonToken {
    enum class Kind : unsigned char { Scalar, Object, Array };
    Kind kind = Kind::Scalar;
    std::string scalar;
    std::vector<std::pair<std::string, JsonToken>> members;
    std::vector<JsonToken> elements;

    JsonToken& field(std::string_view key) {
        for (auto& [name, value] : members)
            if (name == key)
                return value;
        throw std::runtime_error("Missing fixture field: " + std::string{key});
    }
    [[nodiscard]] const JsonToken& field(std::string_view key) const {
        for (const auto& [name, value] : members)
            if (name == key)
                return value;
        throw std::runtime_error("Missing fixture field: " + std::string{key});
    }
    [[nodiscard]] bool contains(std::string_view key) const {
        return std::ranges::any_of(members,
                                   [&](const auto& member) { return member.first == key; });
    }
    void set(std::string key, JsonToken value) {
        for (auto& [name, old] : members)
            if (name == key) {
                old = std::move(value);
                return;
            }
        members.emplace_back(std::move(key), std::move(value));
    }
    void erase(std::string_view key) {
        std::erase_if(members, [&](const auto& member) { return member.first == key; });
    }
    [[nodiscard]] std::string text() const {
        return nlohmann::json::parse(scalar).get<std::string>();
    }
    [[nodiscard]] std::string wire() const {
        if (kind == Kind::Scalar)
            return scalar;
        std::string result = kind == Kind::Object ? "{" : "[";
        bool first = true;
        if (kind == Kind::Object) {
            for (const auto& [key, value] : members) {
                if (!first)
                    result += ',';
                first = false;
                result += nlohmann::json(key).dump() + ":" + value.wire();
            }
        } else {
            for (const auto& value : elements) {
                if (!first)
                    result += ',';
                first = false;
                result += value.wire();
            }
        }
        result += kind == Kind::Object ? '}' : ']';
        return result;
    }
};

class TokenReader {
public:
    explicit TokenReader(std::string_view input) : source_(input) {}
    JsonToken read() {
        skip();
        if (position_ >= source_.size())
            throw std::runtime_error("Incomplete fixture JSON");
        if (source_[position_] == '{') {
            ++position_;
            JsonToken value;
            value.kind = JsonToken::Kind::Object;
            skip();
            while (peek() != '}') {
                const auto key = quoted();
                skip();
                consume(':');
                value.members.emplace_back(nlohmann::json::parse(key).get<std::string>(), read());
                skip();
                if (peek() == '}')
                    break;
                consume(',');
            }
            consume('}');
            return value;
        }
        if (source_[position_] == '[') {
            ++position_;
            JsonToken value;
            value.kind = JsonToken::Kind::Array;
            skip();
            while (peek() != ']') {
                value.elements.push_back(read());
                skip();
                if (peek() == ']')
                    break;
                consume(',');
            }
            consume(']');
            return value;
        }
        JsonToken value;
        if (source_[position_] == '"') {
            value.scalar = quoted();
        } else {
            const auto begin = position_;
            while (position_ < source_.size() && source_[position_] != ',' &&
                   source_[position_] != ']' && source_[position_] != '}' &&
                   source_[position_] != ' ' && source_[position_] != '\n' &&
                   source_[position_] != '\r' && source_[position_] != '\t')
                ++position_;
            if (begin == position_)
                throw std::runtime_error("Empty fixture token");
            value.scalar = source_.substr(begin, position_ - begin);
        }
        return value;
    }

private:
    char peek() {
        skip();
        if (position_ >= source_.size())
            throw std::runtime_error("Missing fixture delimiter");
        return source_[position_];
    }
    void skip() {
        while (position_ < source_.size() &&
               (source_[position_] == ' ' || source_[position_] == '\n' ||
                source_[position_] == '\r' || source_[position_] == '\t'))
            ++position_;
    }
    void consume(char token) {
        if (peek() != token)
            throw std::runtime_error("Unexpected fixture delimiter");
        ++position_;
    }
    std::string quoted() {
        skip();
        const auto begin = position_;
        consume('"');
        while (position_ < source_.size()) {
            const auto token = source_[position_++];
            if (token == '\\') {
                if (position_ == source_.size())
                    throw std::runtime_error("Incomplete fixture escape");
                ++position_;
            } else if (token == '"') {
                return std::string{source_.substr(begin, position_ - begin)};
            }
        }
        throw std::runtime_error("Unterminated fixture string");
    }
    std::string_view source_;
    std::size_t position_ = 0;
};

inline JsonToken parse(std::string_view source) {
    return TokenReader{source}.read();
}
inline JsonToken string(std::string_view text) {
    return parse(nlohmann::json(text).dump());
}
inline JsonToken array() {
    return parse("[]");
}
inline JsonToken fragment(std::string_view name) {
    for (const auto& example : case_schema_corpus::fragments)
        if (example.name == name)
            return parse(example.json);
    throw std::runtime_error("Unknown named fragment");
}
inline JsonToken root(std::string_view name = "structural_static") {
    for (const auto& example : case_schema_corpus::roots)
        if (example.name == name)
            return parse(example.json);
    throw std::runtime_error("Unknown named root");
}
inline JsonToken scale(std::string_view field) {
    auto value = parse(
        R"({"kind":"uniform_components","field":"shell_position","targets":[],"values":[1,1,1]})");
    value.set("field", string(field));
    return value;
}

inline JsonToken physical_case(std::string_view context_wire, std::string_view perturbation_wire) {
    const auto context = parse(context_wire);
    auto value = root();
    value.set("units", context.field("units"));
    // Append the optional root member before borrowing any vector-backed member.
    value.set("initial",
              parse(R"({"kind":"equilibrium","branch":"1","perturbation":{"kind":"none"}})"));
    value.field("initial").set("perturbation", parse(perturbation_wire));
    auto& model = value.field("model");
    auto& analysis = value.field("analysis");
    analysis.set("kind", string("transient"));
    analysis.erase("speed");
    analysis.set("intervals", parse(R"([{"start":0,"end":1,"step":1,"restart":"initial"}])"));
    auto& scales = value.field("numerics").field("scales");
    if (!context.field("shell_ids").elements.empty()) {
        model.field("shells").elements[0].set("id", context.field("shell_ids").elements[0]);
        if (context.contains("reference")) {
            const auto& reference = context.field("reference");
            auto& geometry = model.field("geometry").elements[0];
            geometry.set("n", reference.field("n"));
            geometry.set("x0", reference.field("x0"));
            geometry.set("width", reference.field("width"));
            geometry.set("skew", reference.field("skew"));
            geometry.set("reference_length", reference.field("reference_length"));
            for (auto& thickness : geometry.field("thickness").field("values").elements)
                thickness = reference.field("thickness");
        }
        scales.elements.push_back(scale("structural_velocity"));
    } else {
        model.set("geometry", array());
        model.set("materials", array());
        model.set("shells", array());
        auto beam = fragment("Beam");
        beam.set("id", context.field("beam_ids").elements[0]);
        beam.set("nodes", context.field("beam_reference_frames"));
        beam.field("coordinate").set("length", context.field("beam_length"));
        model.set("beams", array());
        model.field("beams").elements.push_back(beam);
        auto rigid = fragment("Mass.rigid");
        rigid.set("id", context.field("rigid_ids").elements[0]);
        rigid.field("target").field("target").set("kind", string("beam"));
        rigid.field("target").field("target").set("id", context.field("beam_ids").elements[0]);
        rigid.field("target").set("coordinates", parse("[1,0,0]"));
        rigid.field("frame").set("origin", context.field("equilibrium_rigid").field("position"));
        rigid.field("frame").set("rotation", context.field("equilibrium_rigid").field("rotation"));
        rigid.field("ledger").field("owner").set("kind", string("beam"));
        rigid.field("ledger").field("owner").set("id", context.field("beam_ids").elements[0]);
        model.set("masses", array());
        model.field("masses").elements.push_back(std::move(rigid));
        auto support = fragment("Support.constraint");
        support.field("target").field("target").set("kind", string("beam"));
        support.field("target").field("target").set("id", context.field("beam_ids").elements[0]);
        support.field("target").set("coordinates", parse("[0,0,0]"));
        support.set("components", parse("[0,1,2,3,4,5]"));
        model.set("supports", array());
        model.field("supports").elements.push_back(std::move(support));
        value.set("observations", array());
        value.field("output").field("observations").set("observations", array());
        value.field("numerics").set("patch", parse(R"({"kind":"none"})"));
        scales = array();
        for (auto name : {"beam_translation", "beam_rotation", "rigid_translation",
                          "rigid_rotation", "force_balance", "moment_balance",
                          "structural_velocity", "structural_angular_velocity"})
            scales.elements.push_back(scale(name));
        scales.elements.push_back(parse(
            R"({"kind":"physical_layout_policy","field":"mixed_blocks","constraint":{"kind":"physical","position":{"kind":"values","values":[1,1,1]},"orientation":{"kind":"values","values":[1,1,1]},"multiplier_work_scale":1,"reaction_force":[1,1,1],"reaction_moment":[1,1,1]}})"));
    }
    return value;
}
} // namespace uvlas::test::schema_adapter
