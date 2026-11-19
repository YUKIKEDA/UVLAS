#include "reference_declarations.hpp"

#include "case_schema_storage.hpp"

#include <algorithm>
#include <charconv>
#include <limits>
#include <new>
#include <utility>

namespace uvlas::io::detail {
namespace {
constexpr auto operation = core::Operation::PrepareModel;
constexpr std::array<std::string_view, 8> role_names{
    "nastran-deck",       "nastran-sidecar",     "gmsh-msh4.1-ascii", "gmsh-sidecar",
    "reference-geometry", "preparation-history", "initial-history",   "observation-reference"};
core::Failure invalid() noexcept {
    return {core::FailureCode::InvalidInput, operation};
}
struct Range {
    std::uint64_t offset = 0, count = 0;
};
struct Row {
    SourceId id;
    ReferenceRole role;
    SourceDigest digest;
    Range spelling, path, hex, pointer;
    std::uint64_t node;
};
struct Stop {
    core::Failure failure;
};
} // namespace
struct ReferenceDeclarationData {
    std::optional<SchemaIdentity> identity;
    core::PlannedBuffer<std::optional<Row>> rows;
    core::PlannedBuffer<char> text;
    DeclarationMemoryReport report{};
};
struct ReferenceDeclarationFactory {
    static const SchemaData* storage(const SchemaCase& source) noexcept {
        return source.data_.size() ? &source.data_.values()[0] : nullptr;
    }
    static DeclaredReferences result(core::PlannedBuffer<ReferenceDeclarationData> data) noexcept {
        return DeclaredReferences{std::move(data)};
    }
    static void diagnostic(DeclarationDiagnostic& slot, DeclarationDiagnosticDetail detail,
                           core::PlannedBuffer<char> path) noexcept {
        slot.detail_ = detail;
        slot.path_ = std::move(path);
    }
};
namespace {
struct Context {
    const SchemaData* source;
    const DeclarationLimits& limits;
    core::ResourcePlan& plan;
    std::stop_token stop;
    std::optional<SchemaIdentity> identity;
    std::uint64_t start, work = 0;
    std::optional<std::uint64_t> current{};
    DeclarationReason reason = DeclarationReason::UnexpectedException;
    std::optional<DeclarationQuota> quota{};
    std::optional<std::uint64_t> required{}, limit{};

    [[noreturn]] void fail(DeclarationReason r, core::FailureCode code) {
        reason = r;
        throw Stop{{code, operation}};
    }
    void poll() {
        if (stop.stop_requested())
            fail(DeclarationReason::Cancelled, core::FailureCode::Cancelled);
    }
    [[nodiscard]] std::uint64_t cap(DeclarationQuota q) const {
        return *limits.get(q);
    }
    void check(DeclarationQuota q, std::uint64_t n) {
        if (n > cap(q)) {
            quota = q;
            required = n;
            limit = cap(q);
            fail(DeclarationReason::QuotaExceeded, core::FailureCode::ResourceExceeded);
        }
    }
    std::uint64_t add(std::uint64_t a, std::uint64_t b) {
        if (b > UINT64_MAX - a)
            fail(DeclarationReason::Overflow, core::FailureCode::ResourceExceeded);
        return a + b;
    }
    std::uint64_t mul(std::uint64_t a, std::uint64_t b) {
        if (a && b > UINT64_MAX / a)
            fail(DeclarationReason::Overflow, core::FailureCode::ResourceExceeded);
        return a * b;
    }
    void step(std::uint64_t n = 1) {
        poll();
        const auto w = add(work, n);
        check(DeclarationQuota::Work, w);
        work = w;
    }
    [[nodiscard]] const InputDocumentData& input() const {
        return source->snapshot.values()[0];
    }
    [[nodiscard]] const SchemaNode& node(std::uint64_t i) const {
        return source->nodes.values()[static_cast<std::size_t>(i)];
    }
    [[nodiscard]] std::string_view text(TextRange r) const {
        return {input().text.values().data() + r.offset, r.size};
    }
    std::optional<std::uint64_t> find(std::uint64_t object, std::string_view key) {
        const auto& value =
            std::get<ObjectValue>(input().nodes.values()[static_cast<std::size_t>(object)].value);
        for (auto i = value.members.head; i >= 0;) {
            step();
            const auto& member = input().members.values()[static_cast<std::size_t>(i)];
            const auto name = text(member.key);
            step(add(name.size(), key.size()));
            if (name == key)
                return static_cast<std::uint64_t>(member.child);
            i = member.next;
        }
        return std::nullopt;
    }
    std::uint64_t field(std::uint64_t object, std::string_view key) {
        const auto v = find(object, key);
        if (!v)
            fail(DeclarationReason::InvalidOwner, core::FailureCode::InvalidInput);
        return *v;
    }
    [[nodiscard]] std::string_view string(std::uint64_t i) const {
        return node_text(input(), i);
    }
    std::optional<std::uint64_t> child(std::uint64_t parent,
                                       std::optional<std::uint64_t> after = {}) {
        const auto kind = node_kind(input(), parent);
        if (kind == InputNodeKind::Array) {
            const auto index = after ? add(*node(*after).parent_array_index, 1) : 0;
            if (index < array_size(input(), parent)) {
                step();
                return source->array_elements
                    .values()[static_cast<std::size_t>(add(node(parent).array_offset, index))];
            }
        } else if (kind == InputNodeKind::Object) {
            const auto& fields = schema_records()[node(parent).record_kind].fields;
            bool ready = !after;
            for (const auto& f : fields) {
                step();
                if (!ready) {
                    const auto member = *node(*after).parent_member;
                    const auto name =
                        text(input().members.values()[static_cast<std::size_t>(member)].key);
                    step(add(name.size(), f.name.size()));
                    if (name == f.name)
                        ready = true;
                    continue;
                }
                if (auto v = find(parent, f.name))
                    return v;
            }
        }
        return std::nullopt;
    }
    std::uint64_t pointer(std::uint64_t original, char* output = nullptr,
                          std::uint64_t capacity = 0) {
        std::uint64_t size = 0;
        auto cursor = capacity;
        for (auto i = original; node(i).parent; i = *node(i).parent) {
            step();
            const auto& n = node(i);
            std::array<char, 20> digits{};
            std::string_view part;
            if (n.parent_member)
                part =
                    text(input().members.values()[static_cast<std::size_t>(*n.parent_member)].key);
            else {
                const auto value = n.parent_array_index;
                if (!value)
                    fail(DeclarationReason::InvalidOwner, core::FailureCode::InvalidInput);
                const auto end =
                    std::to_chars(digits.data(), digits.data() + digits.size(), *value);
                part = {digits.data(), static_cast<std::size_t>(end.ptr - digits.data())};
                step(part.size());
            }
            std::uint64_t length = 1;
            for (const char ch : part) {
                step();
                length = add(length, ch == '~' || ch == '/' ? 2 : 1);
            }
            size = add(size, length);
            if (output) {
                if (length > cursor)
                    fail(DeclarationReason::InvalidOwner, core::FailureCode::InvalidInput);
                cursor -= length;
                auto at = cursor;
                output[at++] = '/';
                for (const char ch : part) {
                    step();
                    output[at++] = ch == '~' || ch == '/' ? '~' : ch;
                    if (ch == '~' || ch == '/')
                        output[at++] = ch == '~' ? '0' : '1';
                }
            }
        }
        if (output && cursor != 0)
            fail(DeclarationReason::InvalidOwner, core::FailureCode::InvalidInput);
        return size;
    }
};
struct Totals {
    std::uint64_t rows = 0, text = 0, paths = 0, maximum_pointer = 0, depth = 0;
};
struct Capture {
    std::uint64_t id, node;
    ReferenceRole role;
    SourceDigest digest;
    std::string_view spelling, path, hex;
    std::uint64_t pointer;
};
Capture capture(Context& c, std::uint64_t node) {
    const auto id = c.field(node, "source_id");
    const auto path = c.field(node, "path");
    const auto hex = c.field(node, "sha256");
    const auto role = c.string(c.field(node, "role"));
    std::optional<ReferenceRole> found;
    for (std::size_t i = 0; i < role_names.size(); ++i) {
        c.step(c.add(role.size(), role_names[i].size()));
        if (role == role_names[i]) {
            found = static_cast<ReferenceRole>(i);
            break;
        }
    }
    if (!found)
        c.fail(DeclarationReason::UnsupportedRole, core::FailureCode::Unsupported);
    const auto digest = c.string(hex);
    if (digest.size() != 64)
        c.fail(DeclarationReason::InvalidOwner, core::FailureCode::InvalidInput);
    SourceDigest bytes{};
    auto nibble = [&](char x) -> unsigned {
        c.step();
        if (x >= '0' && x <= '9')
            return static_cast<unsigned>(x - '0');
        if (x >= 'a' && x <= 'f')
            return static_cast<unsigned>(x - 'a' + 10);
        c.fail(DeclarationReason::InvalidOwner, core::FailureCode::InvalidInput);
    };
    for (std::size_t i = 0; i < 32; ++i)
        bytes.bytes[i] =
            static_cast<std::byte>((nibble(digest[2 * i]) << 4U) | nibble(digest[2 * i + 1]));
    return {c.node(id).unsigned_value,
            node,
            *found,
            bytes,
            c.string(id),
            c.string(path),
            digest,
            c.pointer(node)};
}
// Parent metadata makes the traversal bounded without recursive or dynamically allocated frames.
// Object keys are looked up in descriptor order; source JSON order is never an ordering key.
template <class Emit>
Totals walk(Context& c, Emit emit) {
    Totals total;
    std::uint64_t current = 0, depth = 0;
    for (;;) {
        c.current = current;
        c.step();
        const auto kind = node_kind(c.input(), current);
        const bool container = kind == InputNodeKind::Object || kind == InputNodeKind::Array;
        total.depth = std::max(total.depth, c.add(depth, container ? 1 : 0));
        c.check(DeclarationQuota::SourceDepth, total.depth);
        bool declaration = false;
        if (kind == InputNodeKind::Object &&
            schema_records()[c.node(current).record_kind].family == "FileReference") {
            declaration = true;
            const auto value = capture(c, current);
            total.rows = c.add(total.rows, 1);
            c.check(DeclarationQuota::Declarations, total.rows);
            total.paths = c.add(total.paths, value.path.size());
            c.check(DeclarationQuota::SourcePathBytes, total.paths);
            total.text = c.add(total.text, c.add(c.add(value.spelling.size(), value.path.size()),
                                                 c.add(value.hex.size(), value.pointer)));
            c.check(DeclarationQuota::OwnedTextBytes, total.text);
            total.maximum_pointer = std::max(total.maximum_pointer, value.pointer);
            emit(value);
        }
        if (!declaration)
            if (auto first = c.child(current)) {
                current = *first;
                depth = c.add(depth, 1);
                continue;
            }
        for (;;) {
            if (!c.node(current).parent)
                return total;
            const auto parent = *c.node(current).parent;
            if (auto next = c.child(parent, current)) {
                current = *next;
                break;
            }
            current = parent;
            --depth;
            c.step();
        }
    }
}
// Automatic C++ value-payload envelope, excluding ABI frame padding, allocator/RSS and callers.
constexpr std::uint64_t scratch_bytes = sizeof(Context) + 2 * sizeof(Totals) + 2 * sizeof(Capture) +
                                        8 * sizeof(std::uint64_t) +
                                        2 * sizeof(std::array<char, 20>) + sizeof(Row);
std::string_view view(const ReferenceDeclarationData& d, Range r) noexcept {
    if (!r.count)
        return {};
    return {d.text.values().data() + r.offset, static_cast<std::size_t>(r.count)};
}
void diagnose(Context& c, DeclarationDiagnostic* slot) noexcept {
    if (!slot)
        return;
    DeclarationDiagnosticDetail fixed{c.reason,
                                      c.identity,
                                      c.identity ? c.current : std::nullopt,
                                      c.quota,
                                      c.required,
                                      c.limit,
                                      c.work,
                                      SchemaPathState::NotBuilt};
    core::PlannedBuffer<char> path;
    try {
        if (c.identity && c.current) {
            const auto n = c.pointer(*c.current);
            const auto spent = c.plan.allocated_bytes() - c.start;
            if (n <= c.cap(DeclarationQuota::DiagnosticBytes) &&
                spent <= c.cap(DeclarationQuota::CumulativeBytes) &&
                n <= c.cap(DeclarationQuota::CumulativeBytes) - spent &&
                n <= c.plan.remaining_bytes() && n <= static_cast<std::uint64_t>(INT64_MAX)) {
                // Precharge the complete ancestry/copy envelope before path allocation.
                c.check(DeclarationQuota::Work,
                        c.add(c.work, c.add(c.mul(c.work - fixed.work, 2), n)));
                c.step(n);
                auto p = c.plan.allocate<char>(static_cast<core::Index>(n), operation);
                if (p) {
                    path = std::move(*p);
                    if (n)
                        c.pointer(*c.current, path.values().data(), n);
                    fixed.path_state = SchemaPathState::Complete;
                }
            }
        }
    } catch (...) {
        // The primary value failure remains authoritative; never expose a partial pointer.
        path = {};
        fixed.path_state = SchemaPathState::NotBuilt;
    }
    fixed.work = c.work;
    ReferenceDeclarationFactory::diagnostic(*slot, fixed, std::move(path));
}
} // namespace
std::optional<std::string_view> reference_role_name(ReferenceRole role) noexcept {
    const auto i = static_cast<std::size_t>(role);
    return i < role_names.size() ? std::optional{role_names[i]} : std::nullopt;
}
core::Result<std::uint64_t> DeclarationLimits::get(DeclarationQuota q) const noexcept {
    const auto i = static_cast<std::size_t>(q);
    if (i >= caps_.size())
        return std::unexpected(invalid());
    return caps_[i];
}
DeclarationDiagnostic::DeclarationDiagnostic(DeclarationDiagnostic&& other) noexcept {
    *this = std::move(other);
}
DeclarationDiagnostic& DeclarationDiagnostic::operator=(DeclarationDiagnostic&& other) noexcept {
    if (this != &other) {
        detail_ = other.detail_;
        path_ = std::move(other.path_);
        other.reset();
    }
    return *this;
}
std::optional<std::string_view> DeclarationDiagnostic::path() const& noexcept {
    if (!detail_ || detail_->path_state != SchemaPathState::Complete)
        return std::nullopt;
    return path_.size() ? std::string_view{path_.values().data(), path_.size()}
                        : std::string_view{};
}
void DeclarationDiagnostic::reset() noexcept {
    detail_.reset();
    path_ = {};
}
DeclaredReferences::DeclaredReferences(core::PlannedBuffer<ReferenceDeclarationData> data) noexcept
    : data_(std::move(data)) {}
DeclaredReferences::DeclaredReferences(DeclaredReferences&&) noexcept = default;
DeclaredReferences& DeclaredReferences::operator=(DeclaredReferences&&) noexcept = default;
DeclaredReferences::~DeclaredReferences() = default;
core::Result<SchemaIdentity> DeclaredReferences::identity() const noexcept {
    if (!data_.size() || !data_.values()[0].identity)
        return std::unexpected(invalid());
    return *data_.values()[0].identity;
}
core::Result<std::uint64_t> DeclaredReferences::size() const noexcept {
    if (!data_.size())
        return std::unexpected(invalid());
    return data_.values()[0].rows.size();
}
core::Result<DeclaredReferenceView> DeclaredReferences::entry(core::Index index) const& noexcept {
    if (!data_.size() || index < 0)
        return std::unexpected(invalid());
    const auto& d = data_.values()[0];
    if (static_cast<std::uint64_t>(index) >= d.rows.size())
        return std::unexpected(invalid());
    const auto& row = d.rows.values()[static_cast<std::size_t>(index)];
    if (!row)
        return std::unexpected(invalid());
    return DeclaredReferenceView{row->id,
                                 row->role,
                                 row->digest,
                                 view(d, row->spelling),
                                 view(d, row->path),
                                 view(d, row->hex),
                                 row->node,
                                 view(d, row->pointer)};
}
core::Result<DeclarationMemoryReport> DeclaredReferences::memory_report() const noexcept {
    if (!data_.size())
        return std::unexpected(invalid());
    return data_.values()[0].report;
}
core::Result<DeclaredReferences>
collect_case_file_references(const SchemaCase& source, const DeclarationLimits& limits,
                             core::ResourcePlan& plan, std::stop_token stop,
                             DeclarationDiagnostic* slot) noexcept {
    if (slot)
        slot->reset();
    Context c{ReferenceDeclarationFactory::storage(source),
              limits,
              plan,
              std::move(stop),
              std::nullopt,
              plan.allocated_bytes()};
    try {
        c.poll();
        const auto identity = source.identity();
        if (!c.source || !identity)
            c.fail(DeclarationReason::InvalidOwner, core::FailureCode::InvalidInput);
        c.identity = *identity;
        c.check(DeclarationQuota::SourceNodes, c.source->nodes.size());
        c.check(DeclarationQuota::ScratchBytes, scratch_bytes);
        const auto totals = walk(c, [](const Capture&) {});
        const auto first_work = c.work;
        c.current.reset();
        const auto owned =
            c.add(sizeof(ReferenceDeclarationData),
                  c.add(c.mul(totals.rows, sizeof(std::optional<Row>)), totals.text));
        c.check(DeclarationQuota::OwnedBytes, owned);
        // Repeat the metadata walk, write the complete text, and account for array initialization.
        // Pointer writer adds at most its earlier sizing work plus output bytes.
        const auto forecast =
            c.add(c.mul(first_work, 3), c.add(c.mul(totals.text, 2), c.mul(owned, 2)));
        c.check(DeclarationQuota::Work, forecast);
        const auto reserve =
            std::min(totals.maximum_pointer, c.cap(DeclarationQuota::DiagnosticBytes));
        const auto admission = c.add(owned, reserve);
        c.check(DeclarationQuota::CumulativeBytes, admission);
        if (admission > plan.remaining_bytes()) {
            c.quota = DeclarationQuota::CumulativeBytes;
            c.required = admission;
            c.limit = plan.remaining_bytes();
            c.fail(DeclarationReason::QuotaExceeded, core::FailureCode::ResourceExceeded);
        }
        if (totals.rows > static_cast<std::uint64_t>(INT64_MAX) ||
            totals.text > static_cast<std::uint64_t>(INT64_MAX))
            c.fail(DeclarationReason::Overflow, core::FailureCode::ResourceExceeded);
        auto allocate = [&]<class T>(std::uint64_t n) {
            c.poll();
            c.step(c.mul(n, sizeof(T)));
            auto p = plan.allocate<T>(static_cast<core::Index>(n), operation);
            if (!p)
                c.fail(DeclarationReason::Allocation, core::FailureCode::ResourceExceeded);
            c.poll();
            return std::move(*p);
        };
        auto output = allocate.template operator()<ReferenceDeclarationData>(1);
        auto& d = output.values()[0];
        d.rows = allocate.template operator()<std::optional<Row>>(totals.rows);
        d.text = allocate.template operator()<char>(totals.text);
        std::uint64_t row = 0, at = 0;
        auto copy = [&](std::string_view s) {
            c.step(s.size());
            Range r{at, s.size()};
            std::copy(s.begin(), s.end(),
                      d.text.values().begin() + static_cast<std::ptrdiff_t>(at));
            at = c.add(at, s.size());
            return r;
        };
        walk(c, [&](const Capture& v) {
            const auto spelling = copy(v.spelling), path = copy(v.path), hex = copy(v.hex);
            const Range pointer{at, v.pointer};
            if (v.pointer)
                c.pointer(v.node, d.text.values().data() + at, v.pointer);
            at = c.add(at, v.pointer);
            c.step(sizeof(Row));
            d.rows.values()[static_cast<std::size_t>(row++)] =
                Row{SourceId{v.id}, v.role, v.digest, spelling, path, hex, pointer, v.node};
        });
        c.current.reset();
        c.step(sizeof(ReferenceDeclarationData));
        c.poll();
        d.identity = *identity;
        d.report = {owned,    scratch_bytes, plan.allocated_bytes() - c.start,
                    c.work,   owned,         scratch_bytes,
                    forecast, admission};
        return ReferenceDeclarationFactory::result(std::move(output));
    } catch (const Stop& error) {
        diagnose(c, slot);
        return std::unexpected(error.failure);
    } catch (const std::bad_alloc&) {
        c.reason = DeclarationReason::Allocation;
        diagnose(c, slot);
        return std::unexpected(core::Failure{core::FailureCode::ResourceExceeded, operation});
    } catch (...) {
        c.reason = DeclarationReason::UnexpectedException;
        diagnose(c, slot);
        return std::unexpected(invalid());
    }
}
} // namespace uvlas::io::detail
