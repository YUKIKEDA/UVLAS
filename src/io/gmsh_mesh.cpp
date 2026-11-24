#include "gmsh_mesh.hpp"

#include <algorithm>
#include <charconv>
#include <cstring>
#include <functional>
#include <limits>
#include <utility>

namespace uvlas::io::detail {
namespace {

using enum core::FailureCode;
constexpr auto operation = core::Operation::PrepareModel;

struct Abort {
    core::Failure failure;
    std::size_t offset = 0;
};

[[noreturn]] void reject(core::FailureCode code, std::size_t offset) {
    throw Abort{{code, operation}, offset};
}

class Work {
public:
    Work(std::uint64_t limit, const std::stop_token& stop) noexcept : limit_(limit), stop_(stop) {}

    void charge(std::uint64_t units, std::size_t offset) {
        if (stop_.stop_requested())
            reject(Cancelled, offset);
        if (units > limit_ - used_)
            reject(ResourceExceeded, offset);
        used_ += units;
    }

    [[nodiscard]] std::uint64_t used() const noexcept {
        return used_;
    }

private:
    std::uint64_t limit_, used_ = 0;
    const std::stop_token& stop_;
};

struct Token {
    std::string_view text;
    std::size_t offset;
    [[nodiscard]] GmshLexeme range() const noexcept {
        return {offset, text.size()};
    }
};

[[nodiscard]] bool whitespace(char c) noexcept {
    return c == ' ' || c == '\t' || c == '\n' || c == '\r';
}

class Input {
public:
    Input(std::string_view bytes, Work& work) noexcept : bytes_(bytes), work_(work) {}

    [[nodiscard]] std::size_t position() const noexcept {
        return cursor_;
    }

    bool more() {
        while (cursor_ < bytes_.size() && whitespace(bytes_[cursor_]))
            advance();
        return cursor_ != bytes_.size();
    }

    Token word() {
        if (!more())
            reject(InvalidInput, cursor_);
        const auto first = cursor_;
        while (cursor_ < bytes_.size() && !whitespace(bytes_[cursor_]))
            advance();
        // Budget the second scan of a lexeme by its numeric parser or comparison.
        work_.charge(cursor_ - first, first);
        return {bytes_.substr(first, cursor_ - first), first};
    }

    Token quoted() {
        if (!more() || bytes_[cursor_] != '"')
            reject(InvalidInput, cursor_);
        advance();
        const auto first = cursor_;
        while (cursor_ < bytes_.size() && bytes_[cursor_] != '"') {
            const auto c = bytes_[cursor_];
            if (c == '\n' || c == '\r' || c == '\0')
                reject(InvalidInput, cursor_);
            advance();
        }
        if (cursor_ == bytes_.size())
            reject(InvalidInput, cursor_);
        const auto length = cursor_ - first;
        advance();
        if (cursor_ < bytes_.size() && !whitespace(bytes_[cursor_]))
            reject(InvalidInput, cursor_);
        return {bytes_.substr(first, length), first};
    }

    void expect(std::string_view spelling) {
        const auto token = word();
        if (token.text != spelling)
            reject(InvalidInput, token.offset);
    }

private:
    void advance() {
        work_.charge(1, cursor_);
        ++cursor_;
    }

    std::string_view bytes_;
    Work& work_;
    std::size_t cursor_ = 0;
};

std::uint64_t unsigned_number(Token token) {
    std::uint64_t value = 0;
    const auto* first = token.text.data();
    const auto* end = first + token.text.size();
    const auto parsed = std::from_chars(first, end, value);
    if (parsed.ec != std::errc{} || parsed.ptr != end)
        reject(InvalidInput, token.offset);
    return value;
}

std::uint64_t positive_number(Token token) {
    const auto value = unsigned_number(token);
    if (value == 0)
        reject(InvalidInput, token.offset);
    return value;
}

std::int32_t entity_tag(Token token) {
    const auto value = positive_number(token);
    if (value > static_cast<std::uint64_t>(std::numeric_limits<std::int32_t>::max()))
        reject(InvalidInput, token.offset);
    return static_cast<std::int32_t>(value);
}

GmshLink boundary_link(Token token, std::uint8_t dimension) {
    const bool reversed = token.text.starts_with('-');
    if (reversed)
        token.text.remove_prefix(1);
    if (token.text.empty())
        reject(InvalidInput, token.offset);
    return {{dimension, entity_tag(token)}, reversed};
}

GmshLexeme decimal(Token token) {
    constexpr auto is_digit = [](char c) { return c >= '0' && c <= '9'; };
    std::size_t i = token.text.starts_with('+') || token.text.starts_with('-') ? 1 : 0;
    auto consume_digits = [&]() {
        const auto start = i;
        while (i < token.text.size() && is_digit(token.text[i]))
            ++i;
        return i != start;
    };
    bool digits = consume_digits();
    if (i < token.text.size() && token.text[i] == '.') {
        ++i;
        digits = consume_digits() || digits;
    }
    if (!digits)
        reject(InvalidInput, token.offset);
    if (i < token.text.size() && (token.text[i] == 'e' || token.text[i] == 'E')) {
        ++i;
        if (i < token.text.size() && (token.text[i] == '+' || token.text[i] == '-'))
            ++i;
        if (!consume_digits())
            reject(InvalidInput, token.offset);
    }
    if (i != token.text.size())
        reject(InvalidInput, token.offset);
    return token.range();
}

std::size_t count(Token token, std::uint64_t limit) {
    const auto value = unsigned_number(token);
    if (value > limit ||
        value > static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max()) ||
        value > std::numeric_limits<std::size_t>::max())
        reject(ResourceExceeded, token.offset);
    return static_cast<std::size_t>(value);
}

std::size_t block_count(Token token, std::size_t remaining) {
    const auto value = unsigned_number(token);
    if (value > remaining)
        reject(InvalidInput, token.offset);
    return static_cast<std::size_t>(value);
}

std::uint64_t remaining_count(std::uint64_t limit, std::size_t used) noexcept {
    const auto addressable = std::min<std::uint64_t>(
        std::numeric_limits<std::size_t>::max(),
        static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max()));
    return std::min(limit, addressable) - used;
}

std::uint8_t dimension(Token token) {
    const auto value = unsigned_number(token);
    if (value > 3)
        reject(InvalidInput, token.offset);
    return static_cast<std::uint8_t>(value);
}

struct ElementType {
    GmshElementKind kind;
    std::uint8_t dimension;
    std::size_t corners;
};

ElementType element_type(Token token) {
    switch (unsigned_number(token)) {
    case 15:
        return {GmshElementKind::Point, 0, 1};
    case 1:
        return {GmshElementKind::Line, 1, 2};
    case 2:
        return {GmshElementKind::Triangle, 2, 3};
    case 3:
        return {GmshElementKind::Quadrangle, 2, 4};
    default:
        reject(Unsupported, token.offset);
    }
}

struct Counts {
    std::size_t names = 0, entities = 0, links = 0, nodes = 0, elements = 0;
    std::size_t blocks = 0;
    bool entities_declared = false;
    bool operator==(const Counts&) const = default;
};

struct BlockOwner {
    GmshKey key;
    std::size_t tag_offset = 0;
};

struct Records {
    std::span<GmshNames> names;
    std::span<GmshEntity> entities;
    std::span<GmshLink> links;
    std::span<GmshNode> nodes;
    std::span<GmshElement> elements;
    std::span<BlockOwner> blocks;
};

template <class T>
void append(std::span<T> output, std::size_t index, const T& item) {
    if (!output.empty())
        output[index] = item;
}

struct TagRange {
    std::uint64_t minimum = std::numeric_limits<std::uint64_t>::max(), maximum = 0;

    void add(std::uint64_t tag) noexcept {
        minimum = std::min(minimum, tag);
        maximum = std::max(maximum, tag);
    }

    void check(std::uint64_t lo, std::uint64_t hi, std::size_t offset) const {
        if (hi != maximum || lo != (maximum == 0 ? 0 : minimum))
            reject(InvalidInput, offset);
    }
};

class Parser {
public:
    Parser(std::string_view text, Work& work, const GmshLimits& limits, Records records) noexcept
        : input_(text, work), limits_(limits), output_(records) {}

    Counts read() {
        input_.expect("$MeshFormat");
        const auto version = input_.word();
        const auto file_type = input_.word();
        const auto size = input_.word();
        const auto wordsize = unsigned_number(size);
        if (version.text != "4.1" || unsigned_number(file_type) != 0 ||
            (wordsize != 4 && wordsize != 8))
            reject(Unsupported, version.offset);
        input_.expect("$EndMeshFormat");
        while (input_.more()) {
            const auto section = input_.word();
            if (section.text == "$PhysicalNames") {
                names();
            } else if (section.text == "$Entities") {
                result_.entities_declared = true;
                entities();
            } else if (section.text == "$Nodes") {
                has_nodes_ = true;
                nodes();
            } else if (section.text == "$Elements") {
                has_elements_ = true;
                elements();
            } else {
                reject(section.text.starts_with('$') && !section.text.starts_with("$End")
                           ? Unsupported
                           : InvalidInput,
                       section.offset);
            }
        }
        if (!has_nodes_ || !has_elements_)
            reject(InvalidInput, input_.position());
        return result_;
    }

private:
    GmshKey block_owner() {
        const auto dim = dimension(input_.word());
        const auto tag = input_.word();
        const GmshKey owner{dim, entity_tag(tag)};
        append(output_.blocks, result_.blocks++, BlockOwner{owner, tag.offset});
        return owner;
    }

    void names() {
        const auto limit =
            count(input_.word(), remaining_count(limits_.physical_names, result_.names));
        for (std::size_t i = 0; i < limit; ++i) {
            const auto dim = dimension(input_.word());
            const auto tag = entity_tag(input_.word());
            const auto name = input_.quoted();
            append(output_.names, result_.names++, GmshNames{{dim, tag}, name.range()});
        }
        input_.expect("$EndPhysicalNames");
    }

    void entities() {
        std::array<std::size_t, 4> dimensions{};
        auto record = result_.entities;
        for (auto& amount : dimensions) {
            amount = count(input_.word(), remaining_count(limits_.entities, result_.entities));
            result_.entities += amount;
        }
        for (std::size_t dim = 0; dim < dimensions.size(); ++dim) {
            for (std::size_t i = 0; i < dimensions[dim]; ++i)
                entity(static_cast<std::uint8_t>(dim), record++);
        }
        input_.expect("$EndEntities");
    }

    void entity(std::uint8_t dim, std::size_t record) {
        GmshEntity value{};
        value.key = {dim, entity_tag(input_.word())};
        const auto entries = dim == 0 ? 3U : 6U;
        for (std::size_t i = 0; i < entries; ++i)
            value.bounds[i] = decimal(input_.word());
        value.first_physical = result_.links;
        value.physical_count =
            count(input_.word(), remaining_count(limits_.entity_links, result_.links));
        for (std::size_t i = 0; i < value.physical_count; ++i) {
            const auto tag = entity_tag(input_.word());
            append(output_.links, result_.links++, GmshLink{{dim, tag}, false});
        }
        value.first_boundary = result_.links;
        if (dim != 0) {
            value.boundary_count =
                count(input_.word(), remaining_count(limits_.entity_links, result_.links));
            for (std::size_t i = 0; i < value.boundary_count; ++i) {
                const auto link = boundary_link(input_.word(), static_cast<std::uint8_t>(dim - 1));
                append(output_.links, result_.links++, link);
            }
        }
        append(output_.entities, record, value);
    }

    void nodes() {
        const auto blocks = unsigned_number(input_.word());
        auto written = result_.nodes;
        result_.nodes += count(input_.word(), remaining_count(limits_.nodes, result_.nodes));
        const auto lo = unsigned_number(input_.word());
        const auto hi = unsigned_number(input_.word());
        TagRange range;
        for (std::uint64_t block = 0; block < blocks; ++block) {
            const auto owner = block_owner();
            const auto parameter = input_.word();
            const auto parametric = unsigned_number(parameter);
            if (parametric > 1)
                reject(InvalidInput, parameter.offset);
            const auto amount = block_count(input_.word(), result_.nodes - written);
            for (std::size_t i = 0; i < amount; ++i) {
                const auto tag = positive_number(input_.word());
                range.add(tag);
                append(output_.nodes, written + i, GmshNode{tag, owner, {}, {}, parametric != 0});
            }
            for (std::size_t i = 0; i < amount; ++i) {
                std::array<GmshLexeme, 3> coordinates{}, parameters{};
                for (auto& coordinate : coordinates)
                    coordinate = decimal(input_.word());
                if (parametric != 0)
                    for (std::size_t k = 0; k < owner.dimension; ++k)
                        parameters[k] = decimal(input_.word());
                if (!output_.nodes.empty()) {
                    auto& node = output_.nodes[written + i];
                    node.coordinates = coordinates;
                    node.parameters = parameters;
                }
            }
            written += amount;
        }
        if (written != result_.nodes)
            reject(InvalidInput, input_.position());
        range.check(lo, hi, input_.position());
        input_.expect("$EndNodes");
    }

    void elements() {
        const auto blocks = unsigned_number(input_.word());
        auto written = result_.elements;
        result_.elements +=
            count(input_.word(), remaining_count(limits_.elements, result_.elements));
        const auto lo = unsigned_number(input_.word());
        const auto hi = unsigned_number(input_.word());
        TagRange range;
        for (std::uint64_t block = 0; block < blocks; ++block) {
            const auto owner = block_owner();
            const auto type = element_type(input_.word());
            if (type.dimension != owner.dimension)
                reject(InvalidInput, input_.position());
            const auto amount = block_count(input_.word(), result_.elements - written);
            for (std::size_t i = 0; i < amount; ++i) {
                const auto tag = positive_number(input_.word());
                range.add(tag);
                GmshElement value{tag, owner, type.kind, {}};
                for (std::size_t corner = 0; corner < type.corners; ++corner)
                    value.nodes[corner] = positive_number(input_.word());
                append(output_.elements, written++, value);
            }
        }
        if (written != result_.elements)
            reject(InvalidInput, input_.position());
        range.check(lo, hi, input_.position());
        input_.expect("$EndElements");
    }

    Input input_;
    const GmshLimits& limits_;
    Records output_;
    Counts result_;
    bool has_nodes_ = false, has_elements_ = false;
};

template <class T, class Key>
void order_unique(std::span<T> values, Key key, Work& work, std::size_t offset) {
    std::ranges::sort(values, [&](const T& a, const T& b) {
        work.charge(1, offset);
        return std::invoke(key, a) < std::invoke(key, b);
    });
    for (std::size_t i = 1; i < values.size(); ++i) {
        work.charge(1, offset);
        if (std::invoke(key, values[i - 1]) == std::invoke(key, values[i]))
            reject(InvalidInput, offset);
    }
}

template <class T, class Key, class Id>
void find(std::span<T> values, Key key, Id target, Work& work, std::size_t offset) {
    const auto item =
        std::lower_bound(values.begin(), values.end(), target, [&](const T& value, const Id& id) {
            work.charge(1, offset);
            return std::invoke(key, value) < id;
        });
    work.charge(1, offset);
    if (item == values.end() || std::invoke(key, *item) != target)
        reject(InvalidInput, offset);
}

void resolve(Records records, bool entities_declared, Work& work, std::size_t end) {
    order_unique(records.names, &GmshNames::group, work, end);
    order_unique(records.entities, &GmshEntity::key, work, end);
    order_unique(records.nodes, &GmshNode::tag, work, end);
    order_unique(records.elements, &GmshElement::tag, work, end);
    if (entities_declared) {
        for (const auto& entity : records.entities) {
            work.charge(1, end);
            const auto boundaries =
                records.links.subspan(entity.first_boundary, entity.boundary_count);
            for (const auto& boundary : boundaries)
                find(records.entities, &GmshEntity::key, boundary.key, work, end);
        }
        for (const auto& block : records.blocks)
            find(records.entities, &GmshEntity::key, block.key, work, block.tag_offset);
    }
    for (const auto& element : records.elements)
        for (const auto node : element.connectivity())
            find(records.nodes, &GmshNode::tag, node, work, end);
}

template <class T>
core::PlannedBuffer<T> allocate(core::ResourcePlan& plan, std::size_t count) {
    if (count > static_cast<std::size_t>(std::numeric_limits<core::Index>::max()))
        reject(ResourceExceeded, 0);
    auto result = plan.allocate<T>(static_cast<core::Index>(count), operation);
    if (!result)
        throw Abort{result.error(), 0};
    return std::move(*result);
}

template <class T>
void admit(std::size_t count, std::uint64_t& capacity) {
    if (count > capacity / sizeof(T))
        reject(ResourceExceeded, 0);
    capacity -= count * sizeof(T);
}

struct MeshRecords {
    core::PlannedBuffer<GmshNames> names;
    core::PlannedBuffer<GmshEntity> entities;
    core::PlannedBuffer<GmshLink> links;
    core::PlannedBuffer<GmshNode> nodes;
    core::PlannedBuffer<GmshElement> elements;
};

[[nodiscard]] std::size_t block_owners(const Counts& shape) noexcept {
    return shape.entities_declared ? shape.blocks : 0;
}

void admit_records(const Counts& shape, std::uint64_t capacity) {
    admit<GmshNames>(shape.names, capacity);
    admit<GmshEntity>(shape.entities, capacity);
    admit<GmshLink>(shape.links, capacity);
    admit<GmshNode>(shape.nodes, capacity);
    admit<GmshElement>(shape.elements, capacity);
    admit<BlockOwner>(block_owners(shape), capacity);
}

[[nodiscard]] MeshRecords decode(std::string_view input, const Counts& shape, Work& work,
                                 GmshLimits limits, core::ResourcePlan& plan) {
    MeshRecords result{allocate<GmshNames>(plan, shape.names),
                       allocate<GmshEntity>(plan, shape.entities),
                       allocate<GmshLink>(plan, shape.links), allocate<GmshNode>(plan, shape.nodes),
                       allocate<GmshElement>(plan, shape.elements)};
    auto owners = allocate<BlockOwner>(plan, block_owners(shape));
    const Records records{result.names.values(), result.entities.values(), result.links.values(),
                          result.nodes.values(), result.elements.values(), owners.values()};
    if (Parser{input, work, limits, records}.read() != shape)
        reject(InvalidInput, 0);
    resolve(records, shape.entities_declared, work, input.size());
    work.charge(0, input.size());
    return result;
}

[[nodiscard]] ReferenceFileView mesh_source(const AcquiredRootSnapshot& snapshot, GmshLimits limits,
                                            Work& work) {
    work.charge(0, 0);
    const auto source = snapshot.file.view();
    if (!source || source->bytes.empty() || snapshot.location.role != ReferenceRole::GmshMesh)
        reject(InvalidInput, 0);
    if (source->bytes.size() > limits.file_bytes)
        reject(ResourceExceeded, 0);
    return *source;
}

[[nodiscard]] std::string_view text(std::span<const std::byte> bytes) noexcept {
    return {reinterpret_cast<const char*>(bytes.data()), bytes.size()};
}

template <class T, class Parse>
[[nodiscard]] core::Result<T> with_diagnostic(std::size_t* failure_offset, Parse parse) noexcept {
    if (failure_offset != nullptr)
        *failure_offset = 0;
    return core::guard_entry(operation, InvalidInput, {}, [&]() -> core::Result<T> {
        try {
            return parse();
        } catch (const Abort& error) {
            if (failure_offset != nullptr)
                *failure_offset = error.offset;
            return std::unexpected(error.failure);
        }
    });
}
} // namespace

std::span<const std::uint64_t> GmshElement::connectivity() const& noexcept {
    switch (kind) {
    case GmshElementKind::Point:
        return std::span{nodes}.first(1);
    case GmshElementKind::Line:
        return std::span{nodes}.first(2);
    case GmshElementKind::Triangle:
        return std::span{nodes}.first(3);
    case GmshElementKind::Quadrangle:
        return nodes;
    }
    return {};
}

std::string_view GmshMeshView::text(GmshLexeme value) const noexcept {
    if (value.offset > bytes.size() || value.length > bytes.size() - value.offset)
        return {};
    return bytes.substr(value.offset, value.length);
}

core::Result<GmshMeshLeaf> prove_gmsh_mesh_leaf(const AcquiredRootSnapshot& snapshot,
                                                GmshLimits limits, core::ResourcePlan& plan,
                                                const std::stop_token& stop,
                                                std::size_t* failure_offset) noexcept {
    return with_diagnostic<GmshMeshLeaf>(failure_offset, [&] {
        Work work{limits.work, stop};
        const auto source = mesh_source(snapshot, limits, work);
        const auto input = text(source.bytes);
        const auto shape = Parser{input, work, limits, {}}.read();
        admit_records(shape, plan.remaining_bytes());
        work.charge(0, 0);
        // Every supported section consumes only in-file entities, tags or metadata.
        // The parser rejects unknown sections, so a decoded file has no child declarations.
        static_cast<void>(decode(input, shape, work, limits, plan));
        return GmshMeshLeaf{snapshot.location, source.digest, source.bytes.size(), work.used()};
    });
}

core::Result<GmshMesh> GmshMesh::prepare(const AcquiredRootSnapshot& snapshot, GmshLimits limits,
                                         core::ResourcePlan& plan, const std::stop_token& stop,
                                         std::size_t* failure_offset) noexcept {
    return with_diagnostic<GmshMesh>(failure_offset, [&] {
        Work work{limits.work, stop};
        const auto source = mesh_source(snapshot, limits, work);
        const auto input = text(source.bytes);
        const auto shape = Parser{input, work, limits, {}}.read();
        auto capacity = plan.remaining_bytes();
        admit<char>(input.size(), capacity);
        admit_records(shape, capacity);
        work.charge(0, 0);
        GmshMesh result;
        result.location_ = snapshot.location;
        result.digest_ = source.digest;
        result.bytes_ = allocate<char>(plan, input.size());
        const auto copy = result.bytes_.values();
        for (std::size_t offset = 0; offset < input.size();) {
            const auto part = std::min<std::size_t>(65536, input.size() - offset);
            work.charge(part, offset);
            std::memcpy(copy.data() + offset, input.data() + offset, part);
            offset += part;
        }
        auto records = decode({copy.data(), copy.size()}, shape, work, limits, plan);
        result.names_ = std::move(records.names);
        result.entities_ = std::move(records.entities);
        result.links_ = std::move(records.links);
        result.nodes_ = std::move(records.nodes);
        result.elements_ = std::move(records.elements);
        result.entities_declared_ = shape.entities_declared;
        result.work_ = work.used();
        return result;
    });
}

std::optional<GmshMeshView> GmshMesh::view() const& noexcept {
    if (!location_ || bytes_.size() == 0)
        return std::nullopt;
    const auto bytes = bytes_.values();
    return GmshMeshView{*location_,
                        digest_,
                        {bytes.data(), bytes.size()},
                        names_.values(),
                        entities_.values(),
                        links_.values(),
                        nodes_.values(),
                        elements_.values(),
                        entities_declared_,
                        work_};
}

} // namespace uvlas::io::detail
