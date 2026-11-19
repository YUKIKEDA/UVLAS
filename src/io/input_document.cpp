#include "input_document_detail.hpp"
#include "source_digest.hpp"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <bit>
#include <climits>
#include <cmath>
#include <cstring>
#include <limits>
#include <string>
#include <utility>

namespace uvlas::io {
namespace {
constexpr auto operation = core::Operation::PrepareModel;
using detail::InputDocumentData;
using detail::InputDocumentFactory;
using Json = nlohmann::json;
struct Rejection {
    core::Failure failure;
};
struct Context {
    SourceId source;
    std::optional<InputDiagnostic>& diagnostic;
    core::FailureDetails details;
    std::optional<SourceDigest> digest;
    std::optional<ParserMemoryPlan> plan;
    std::optional<std::uint64_t> event;
    std::optional<std::uint64_t> position;

    InputDiagnostic& describe(InputDiagnosticReason reason) noexcept {
        auto& owned = diagnostic.emplace(reason, source);
        owned.digest = digest;
        owned.parser_plan = plan;
        owned.event_ordinal = event;
        owned.byte_offset = position;
        return owned;
    }
    [[noreturn]] void reject(core::FailureCode code, InputDiagnosticReason reason) {
        describe(reason);
        throw Rejection{core::Failure{code, operation, details}};
    }
    [[noreturn]] void quota(InputQuota kind, std::uint64_t limit,
                            std::optional<std::uint64_t> required) {
        auto& owned = describe(InputDiagnosticReason::QuotaExceeded);
        owned.quota = InputQuotaDiagnostic{kind, limit, required};
        throw Rejection{core::Failure{core::FailureCode::ResourceExceeded, operation, details}};
    }
};
std::uint64_t add(std::uint64_t a, std::uint64_t b, Context& context) {
    if (b > std::numeric_limits<std::uint64_t>::max() - a)
        context.reject(core::FailureCode::ResourceExceeded, InputDiagnosticReason::DimensionRange);
    return a + b;
}
std::uint64_t multiply(std::uint64_t a, std::uint64_t b, Context& context) {
    if (a != 0 && b > std::numeric_limits<std::uint64_t>::max() / a)
        context.reject(core::FailureCode::ResourceExceeded, InputDiagnosticReason::DimensionRange);
    return a * b;
}
core::Index as_index(std::uint64_t value, Context& context) {
    if (value > static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max()) ||
        value > std::numeric_limits<std::size_t>::max())
        context.reject(core::FailureCode::ResourceExceeded, InputDiagnosticReason::DimensionRange);
    return static_cast<core::Index>(value);
}
core::Failure complete(core::Failure failure, const core::FailureDetails& known) noexcept {
    auto& inner = failure.details;
    if (!inner.model)
        inner.model = known.model;
    if (!inner.parent)
        inner.parent = known.parent;
    if (!inner.trial)
        inner.trial = known.trial;
    if (!inner.numerical_generation)
        inner.numerical_generation = known.numerical_generation;
    if (!inner.component)
        inner.component = known.component;
    if (!inner.problem)
        inner.problem = known.problem;
    if (!inner.resources)
        inner.resources = known.resources;
    if (!inner.last_accepted)
        inner.last_accepted = known.last_accepted;
    return failure;
}
template <class T>
T take(core::Result<T> value, Context& context) {
    if (!value) {
        context.describe(value.error().code == core::FailureCode::ResourceExceeded
                             ? InputDiagnosticReason::Allocation
                             : InputDiagnosticReason::UnexpectedException);
        throw Rejection{complete(value.error(), context.details)};
    }
    return std::move(*value);
}
consteval SourceDigest fixed_digest(std::string_view text) {
    if (text.size() != 64)
        throw "Invalid fixed source digest";
    SourceDigest digest{};
    const auto digit = [](char c) -> unsigned {
        if (c >= '0' && c <= '9')
            return static_cast<unsigned>(c - '0');
        if (c >= 'a' && c <= 'f')
            return static_cast<unsigned>(c - 'a') + 10U;
        throw "Invalid fixed source digest digit";
    };
    for (std::size_t i = 0; i < 32; ++i)
        digest.bytes[i] = static_cast<std::byte>(16U * digit(text[2 * i]) + digit(text[2 * i + 1]));
    return digest;
}
ParserEnvironment environment(Context& context) {
#if !defined(UVLAS_NLOHMANN_JSON_ALLOCATION_REPAIR) || UVLAS_NLOHMANN_JSON_ALLOCATION_REPAIR != 1
    context.reject(core::FailureCode::Unsupported, InputDiagnosticReason::UnsupportedMemoryModel);
#endif
    if constexpr (CHAR_BIT != 8 || sizeof(void*) != 8 || sizeof(std::size_t) != 8)
        context.reject(core::FailureCode::Unsupported,
                       InputDiagnosticReason::UnsupportedMemoryModel);
    if constexpr (NLOHMANN_JSON_VERSION_MAJOR != 3 || NLOHMANN_JSON_VERSION_MINOR != 12 ||
                  NLOHMANN_JSON_VERSION_PATCH != 0 || JSON_DIAGNOSTICS != 0 ||
                  JSON_DIAGNOSTIC_POSITIONS != 0)
        context.reject(core::FailureCode::Unsupported,
                       InputDiagnosticReason::UnsupportedMemoryModel);
#if defined(_MSVC_STL_VERSION) && defined(_MSVC_STL_UPDATE)
    if constexpr (_MSVC_STL_VERSION != 145 || _MSVC_STL_UPDATE != 202604L ||
                  (_ITERATOR_DEBUG_LEVEL != 0 && _ITERATOR_DEBUG_LEVEL != 2))
        context.reject(core::FailureCode::Unsupported,
                       InputDiagnosticReason::UnsupportedMemoryModel);
    return {ParserLibraryFamily::MsvcStl, _MSVC_STL_UPDATE, 1, 3, 12, 0, 1, 64, 64, 32, 8,
            _ITERATOR_DEBUG_LEVEL};
#elif defined(_GLIBCXX_RELEASE) && defined(__GLIBCXX__) && defined(_GLIBCXX_USE_CXX11_ABI)
#if defined(_GLIBCXX_DEBUG) || defined(_GLIBCXX_PARALLEL) || defined(_GLIBCXX_PROFILE)
    context.reject(core::FailureCode::Unsupported, InputDiagnosticReason::UnsupportedMemoryModel);
#endif
    // Add another revision only after inspecting its source/model correspondence.
    if constexpr (_GLIBCXX_RELEASE != 14 || (__GLIBCXX__ != 20240801 && __GLIBCXX__ != 20240908) ||
                  _GLIBCXX_USE_CXX11_ABI != 1 || sizeof(unsigned long) != 8)
        context.reject(core::FailureCode::Unsupported,
                       InputDiagnosticReason::UnsupportedMemoryModel);
    return {ParserLibraryFamily::Libstdcxx, __GLIBCXX__, 1, 3, 12, 0, 1, 64, 64, 64, 8, 0};
#else
    context.reject(core::FailureCode::Unsupported, InputDiagnosticReason::UnsupportedMemoryModel);
#endif
}
ParserMemoryPlan memory_plan(std::uint64_t n, std::uint64_t depth, Context& context) {
    const auto platform = environment(context);
    const bool msvc = platform.family == ParserLibraryFamily::MsvcStl;
    const std::uint64_t proxy = msvc && platform.iterator_debug_level != 0 ? 16 : 0;
    const auto s = [&](std::uint64_t length) {
        return add(add(multiply(3, length, context), msvc ? 32 : 2, context), proxy, context);
    };
    const auto r = [&](std::uint64_t length) {
        return add(length, msvc ? 1 : 1 + 24 + 4096, context);
    };
    const auto k = multiply(8, n, context);
    const auto count_s = [&](std::uint64_t count, std::uint64_t length) {
        return multiply(count, s(length), context);
    };
    const auto total = [&](std::initializer_list<std::uint64_t> terms) {
        std::uint64_t result = 0;
        for (const auto value : terms)
            result = add(result, value, context);
        return result;
    };
    // Source-derived live slots, including copy-prvalue results of the repair.
    const auto syntax =
        total({count_s(4, k), count_s(2, add(k, 92, context)), count_s(2, add(k, 170, context)),
               count_s(2, 31), count_s(2, 32), s(16), count_s(2, 34), s(12), s(3), count_s(2, 58),
               count_s(2, 20), s(0), count_s(2, add(k, 275, context)), r(add(k, 275, context))});
    const auto range = total({count_s(4, k), count_s(2, add(k, 26, context)), count_s(2, 34), s(12),
                              s(3), s(0), count_s(2, add(k, 60, context)), r(add(k, 60, context))});
    const auto d = std::min(n, depth);
    const auto word = static_cast<std::uint64_t>(platform.word_bits);
    const auto words = add(d / word, d % word == 0 ? 0 : 1, context);
    const ParserMemoryBreakdown breakdown{
        add(multiply(3, n, context), proxy, context), s(n),
        add(multiply(multiply(3, words, context), word / 8, context), multiply(2, proxy, context),
            context),
        syntax, range};
    const auto estimate = total({breakdown.token_vector_bytes, breakdown.token_buffer_bytes,
                                 breakdown.state_stack_bytes, std::max(syntax, range)});
    std::array<ParserSourceEvidence, 12> sources{};
    std::size_t used = 0;
    const auto evidence = [&](ParserSourceRole role, SourceDigest sha) {
        sources[used++] = {role, sha};
    };
    evidence(ParserSourceRole::Parser,
             fixed_digest("8f219283e1b1a5ad3dc413acfa33c4ee3901c0366f60676c3a93b0d37aeda7b6"));
    evidence(ParserSourceRole::Lexer,
             fixed_digest("9ab10d7a6060e048842465836c9f11ec3b064125197f388e2e9792f7ca1211f2"));
    evidence(ParserSourceRole::Exceptions,
             fixed_digest("5c9b21e0a13ac0bab8dea676b9804d739dfd2adc24f036f97ab243ef70e2aae7"));
    evidence(ParserSourceRole::StringConcat,
             fixed_digest("666b4e48add7381f4b995cd451d3923c1a6b261d633dd9c9ad6f80e54298020b"));
    if (msvc) {
        evidence(ParserSourceRole::Vector,
                 fixed_digest("b254d6fc250164a9f4e71ea66e2c1e9401cb2bf8b0056550fd92135518450b51"));
        evidence(ParserSourceRole::String,
                 fixed_digest("ab80a6b60196246d1c8653e6d55f78e3cbaaf4207860a37cf6996cdc833eea7a"));
        evidence(ParserSourceRole::ExceptionRuntime,
                 fixed_digest("cc401650d7cd6fdefbfaec65b6cb83e504b2617a89bae3a42871ac9960eb09c3"));
        evidence(ParserSourceRole::ExceptionHeader,
                 fixed_digest("0dd3dcc67b655fbfc96708f131c9e7f86416d2db95fd33b1e7acc8489f2c32d4"));
        evidence(ParserSourceRole::IteratorStorage,
                 fixed_digest("0409c9668d5bae4abeeceeef3ee29353c1f4039dc5a78565acfe1cb34d442f38"));
    } else {
        evidence(ParserSourceRole::Vector,
                 fixed_digest("b896547442cc7fa6e6eea35f94a64e0b32431e1e3d8abee301bd5aa1f5ff8c0f"));
        evidence(ParserSourceRole::VectorGrowth,
                 fixed_digest("dc5dad860d9f9add7ccbbea1e67f315c6f0f74d6398ccdd6459391a976cd0dc2"));
        evidence(ParserSourceRole::BitVector,
#if defined(__GLIBCXX__) && __GLIBCXX__ == 20240908
                 fixed_digest("7abd919ef3579f6199fb87981f80239e25ec57bc0f80e16e0ca792d44bb14f29"));
#else
                 fixed_digest("a5e4a8743d833591f6809d5d2bdb9998f7f681c1248a5e2c49623676d29ab4ea"));
#endif
        evidence(ParserSourceRole::String,
                 fixed_digest("4f8dbb49a708622d994acc0b6ad3e69cf72a746732c4ff2362d078b4699d478d"));
        evidence(ParserSourceRole::CowString,
                 fixed_digest("0cb671deeb0d99a38b6171290d056034ef6c2c804078a4fca30b4541f16f4d5b"));
        evidence(ParserSourceRole::ExceptionRuntime,
                 fixed_digest("c40ec114b10ec3a67067d2e5f206cb31d4441afb713bb4d8789d7ec85e107f65"));
    }
    return InputDocumentFactory::plan(platform, breakdown, estimate, sources, used);
}
std::string_view text(const InputDocumentData& data, detail::TextRange value) noexcept {
    return {data.text.values().data() + value.offset, value.size};
}
std::string_view member_key(const InputDocumentData& data, core::Index member) noexcept {
    return text(data, data.members.values()[static_cast<std::size_t>(member)].key);
}
const detail::ChildList& array_list(const InputDocumentData& data, core::Index node) noexcept {
    // Only the private factory can create this view, after checking its node kind.
    return std::get_if<detail::ArrayValue>(
               &data.nodes.values()[static_cast<std::size_t>(node)].value)
        ->children;
}
const detail::ChildList& object_list(const InputDocumentData& data, core::Index node) noexcept {
    return std::get_if<detail::ObjectValue>(
               &data.nodes.values()[static_cast<std::size_t>(node)].value)
        ->members;
}
} // namespace

InputLimits::InputLimits(std::uint64_t source_bytes, std::uint64_t depth, std::uint64_t nodes,
                         std::uint64_t members, std::uint64_t decoded_string_bytes,
                         std::uint64_t decimal_token_bytes, std::uint64_t total_text_bytes,
                         ParserAllocationPolicy parser_memory) noexcept
    : maximum_source_bytes(source_bytes),
      maximum_depth(depth),
      maximum_nodes(nodes),
      maximum_members(members),
      maximum_decoded_string_bytes(decoded_string_bytes),
      maximum_decimal_token_bytes(decimal_token_bytes),
      maximum_total_text_bytes(total_text_bytes),
      external_memory(parser_memory) {}
InputDiagnostic::InputDiagnostic(InputDiagnosticReason why, SourceId origin) noexcept
    : reason(why), source(origin) {}
ParserMemoryPlan::ParserMemoryPlan(ParserEnvironment platform, ParserMemoryBreakdown memory,
                                   std::uint64_t estimate,
                                   std::array<ParserSourceEvidence, 12> sources,
                                   std::size_t count) noexcept
    : environment_(platform),
      breakdown_(memory),
      estimate_(estimate),
      sources_(sources),
      source_count_(count) {}
std::uint64_t ParserMemoryPlan::estimate_bytes() const noexcept {
    return estimate_;
}
std::optional<std::uint64_t> ParserMemoryPlan::proven_upper_bytes() const noexcept {
    return std::nullopt;
}
ParserEnvironment ParserMemoryPlan::environment() const noexcept {
    return environment_;
}
ParserMemoryBreakdown ParserMemoryPlan::breakdown() const noexcept {
    return breakdown_;
}
std::span<const ParserSourceEvidence> ParserMemoryPlan::sources() const& noexcept {
    return {sources_.data(), source_count_};
}
namespace detail {
ParserFailureClassification classify_parser_failure(ParserFailureClass kind, int id,
                                                    std::string_view message) noexcept {
    if (kind != ParserFailureClass::Syntax && kind != ParserFailureClass::NumberRange)
        return {core::FailureCode::InvalidInput, InputDiagnosticReason::UnexpectedException};
    if (message == "Unknown exception")
        return {core::FailureCode::ResourceExceeded, InputDiagnosticReason::Allocation};
    if (!message.starts_with("[json.exception."))
        return {core::FailureCode::InvalidInput, InputDiagnosticReason::UnexpectedException};
    if (kind == ParserFailureClass::NumberRange && id == 406)
        return {core::FailureCode::Nonfinite, InputDiagnosticReason::NumberRange};
    if (kind != ParserFailureClass::Syntax)
        return {core::FailureCode::InvalidInput, InputDiagnosticReason::UnexpectedException};
    // The pinned parser appends this delimiter before quoting user-controlled
    // last-read bytes. Only its preceding generated diagnostic is a cause;
    // input text may itself contain delimiters, quotes and diagnostic words.
    const auto cause = message.substr(0, message.find("; last read: '"));
    if (cause.find("ill-formed UTF-8") != std::string_view::npos)
        return {core::FailureCode::InvalidInput, InputDiagnosticReason::InvalidUtf8};
    if (cause.find("surrogate") != std::string_view::npos ||
        cause.find("4 hex digits") != std::string_view::npos ||
        cause.find("after backslash") != std::string_view::npos)
        return {core::FailureCode::InvalidInput, InputDiagnosticReason::InvalidStringEscape};
    return {core::FailureCode::InvalidInput, InputDiagnosticReason::JsonSyntax};
}
InputNodeOrdinal InputDocumentFactory::ordinal(core::Index value) noexcept {
    return InputNodeOrdinal{value};
}
InputNodeView InputDocumentFactory::node(const InputDocumentData* data,
                                         core::Index index) noexcept {
    return InputNodeView{data, index};
}
InputArrayView InputDocumentFactory::array(const InputDocumentData* data,
                                           core::Index index) noexcept {
    return InputArrayView{data, index};
}
InputObjectView InputDocumentFactory::object(const InputDocumentData* data,
                                             core::Index index) noexcept {
    return InputObjectView{data, index};
}
InputDocument InputDocumentFactory::document(core::PlannedBuffer<InputDocumentData> data) noexcept {
    return InputDocument{std::move(data)};
}
ParserMemoryPlan InputDocumentFactory::plan(ParserEnvironment platform,
                                            ParserMemoryBreakdown memory, std::uint64_t estimate,
                                            std::array<ParserSourceEvidence, 12> sources,
                                            std::size_t count) noexcept {
    return ParserMemoryPlan{platform, memory, estimate, sources, count};
}
} // namespace detail
InputNodeView::InputNodeView(const InputDocumentData* data, core::Index index) noexcept
    : owner_(data), node_(index) {}
InputArrayView::InputArrayView(const InputDocumentData* data, core::Index index) noexcept
    : owner_(data), node_(index) {}
InputObjectView::InputObjectView(const InputDocumentData* data, core::Index index) noexcept
    : owner_(data), node_(index) {}
InputNodeOrdinal InputNodeView::ordinal() const noexcept {
    return InputDocumentFactory::ordinal(node_);
}
InputNodeKind InputNodeView::kind() const noexcept {
    return static_cast<InputNodeKind>(
        owner_->nodes.values()[static_cast<std::size_t>(node_)].value.index());
}
std::optional<bool> InputNodeView::boolean() const noexcept {
    const auto* value =
        std::get_if<bool>(&owner_->nodes.values()[static_cast<std::size_t>(node_)].value);
    return value ? std::optional{*value} : std::nullopt;
}
std::optional<std::int64_t> InputNodeView::signed_integer() const noexcept {
    const auto* value =
        std::get_if<std::int64_t>(&owner_->nodes.values()[static_cast<std::size_t>(node_)].value);
    return value ? std::optional{*value} : std::nullopt;
}
std::optional<std::uint64_t> InputNodeView::unsigned_integer() const noexcept {
    const auto* value =
        std::get_if<std::uint64_t>(&owner_->nodes.values()[static_cast<std::size_t>(node_)].value);
    return value ? std::optional{*value} : std::nullopt;
}
std::optional<NativeInteger> InputNodeView::native_integer() const noexcept {
    if (const auto value = signed_integer())
        return NativeInteger{*value};
    if (const auto value = unsigned_integer())
        return NativeInteger{*value};
    return std::nullopt;
}
std::optional<std::string_view> InputNodeView::string() const& noexcept {
    const auto* value = std::get_if<detail::StringValue>(
        &owner_->nodes.values()[static_cast<std::size_t>(node_)].value);
    return value ? std::optional{text(*owner_, value->text)} : std::nullopt;
}
std::optional<std::string_view> InputNodeView::decimal_token() const& noexcept {
    const auto* value = std::get_if<detail::DecimalValue>(
        &owner_->nodes.values()[static_cast<std::size_t>(node_)].value);
    return value ? std::optional{text(*owner_, value->text)} : std::nullopt;
}
std::optional<InputArrayView> InputNodeView::array() const& noexcept {
    return kind() == InputNodeKind::Array
               ? std::optional{InputDocumentFactory::array(owner_, node_)}
               : std::nullopt;
}
std::optional<InputObjectView> InputNodeView::object() const& noexcept {
    return kind() == InputNodeKind::Object
               ? std::optional{InputDocumentFactory::object(owner_, node_)}
               : std::nullopt;
}
core::Index InputArrayView::size() const noexcept {
    return array_list(*owner_, node_).size;
}
core::Index InputObjectView::size() const noexcept {
    return object_list(*owner_, node_).size;
}
core::Result<InputNodeView> InputArrayView::at(core::Index index) const& noexcept {
    const auto& list = array_list(*owner_, node_);
    if (index < 0 || index >= list.size)
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, operation});
    auto edge = list.head;
    for (core::Index i = 0; i < index; ++i)
        edge = owner_->edges.values()[static_cast<std::size_t>(edge)].next;
    return InputDocumentFactory::node(owner_,
                                      owner_->edges.values()[static_cast<std::size_t>(edge)].child);
}
core::Result<InputMemberView> InputObjectView::at(core::Index index) const& noexcept {
    const auto& list = object_list(*owner_, node_);
    if (index < 0 || index >= list.size)
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, operation});
    auto member = list.head;
    for (core::Index i = 0; i < index; ++i)
        member = owner_->members.values()[static_cast<std::size_t>(member)].next;
    return InputMemberView{
        member_key(*owner_, member),
        InputDocumentFactory::node(
            owner_, owner_->members.values()[static_cast<std::size_t>(member)].child)};
}
std::optional<InputNodeView> InputObjectView::find(std::string_view key) const& noexcept {
    auto member = object_list(*owner_, node_).head;
    while (member != -1) {
        const auto& value = owner_->members.values()[static_cast<std::size_t>(member)];
        if (member_key(*owner_, member) == key)
            return InputDocumentFactory::node(owner_, value.child);
        member = value.next;
    }
    return std::nullopt;
}
InputDocument::InputDocument(core::PlannedBuffer<InputDocumentData> data) noexcept
    : data_(std::move(data)) {}
InputDocument::InputDocument(InputDocument&&) noexcept = default;
InputDocument& InputDocument::operator=(InputDocument&&) noexcept = default;
InputDocument::~InputDocument() = default;
core::Result<InputDocumentView> InputDocument::view() const& noexcept {
    if (data_.size() != 1)
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, operation});
    const auto& data = data_.values()[0];
    if (!data.identity || !data.parser_plan || data.node_count == 0)
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, operation});
    return InputDocumentView{*data.identity, std::as_bytes(data.source.values()),
                             InputDocumentFactory::node(&data, 0), data.policy, *data.parser_plan};
}

namespace {
class DocumentSax final : public nlohmann::json_sax<Json> {
public:
    DocumentSax(InputDocumentData& data, core::PlannedBuffer<detail::BuildFrame>& frames,
                const InputLimits& limits, Context& context) noexcept
        : data_(data), frames_(frames), limits_(limits), context_(context) {}

    bool null() override {
        event();
        (void)node(std::monostate{});
        return true;
    }
    bool boolean(bool value) override {
        event();
        (void)node(value);
        return true;
    }
    bool number_integer(Json::number_integer_t value) override {
        event();
        (void)node(value);
        return true;
    }
    bool number_unsigned(Json::number_unsigned_t value) override {
        event();
        (void)node(value);
        return true;
    }
    bool number_float(Json::number_float_t value, const Json::string_t& raw) override {
        event();
        if (!std::isfinite(value))
            context_.reject(core::FailureCode::Nonfinite, InputDiagnosticReason::NumberRange);
        ensure_node();
        (void)node(detail::DecimalValue{store(raw, true)});
        return true;
    }
    bool string(Json::string_t& value) override {
        event();
        ensure_node();
        (void)node(detail::StringValue{store(value, false)});
        return true;
    }
    bool binary(Json::binary_t&) override {
        event();
        context_.reject(core::FailureCode::InvalidInput,
                        InputDiagnosticReason::UnexpectedException);
    }
    bool start_object(std::size_t) override {
        event();
        start(detail::ObjectValue{});
        return true;
    }
    bool start_array(std::size_t) override {
        event();
        start(detail::ArrayValue{});
        return true;
    }
    bool end_object() override {
        event();
        end(false);
        return true;
    }
    bool end_array() override {
        event();
        end(true);
        return true;
    }
    bool key(Json::string_t& key_value) override {
        event();
        if (depth_ == 0)
            impossible();
        auto& frame = frames_.values()[depth_ - 1];
        auto* object = std::get_if<detail::ObjectValue>(&at(frame.node).value);
        if (object == nullptr || frame.pending_member != -1)
            impossible();
        const auto required = add(data_.member_count, 1, context_);
        if (required > limits_.maximum_members)
            context_.quota(InputQuota::Members, limits_.maximum_members, required);
        check_text(key_value, false);
        const auto hash = key_hash(frame.node, key_value);
        auto slot = static_cast<std::size_t>(hash % data_.keys.size());
        for (std::size_t inspected = 0; inspected < data_.keys.size(); ++inspected) {
            const auto& existing = data_.keys.values()[slot];
            if (existing.member == -1)
                break;
            if (existing.parent == frame.node && existing.hash == hash &&
                member_key(data_, existing.member) == key_value)
                context_.reject(core::FailureCode::InvalidInput,
                                InputDiagnosticReason::DuplicateKey);
            slot = slot + 1 == data_.keys.size() ? 0 : slot + 1;
        }
        if (data_.keys.values()[slot].member != -1 || data_.member_count >= data_.members.size())
            impossible();
        const auto index = static_cast<core::Index>(data_.member_count++);
        auto& member = data_.members.values()[static_cast<std::size_t>(index)];
        member.key = store(key_value, false);
        member.event = current_event();
        data_.keys.values()[slot] = {frame.node, index, hash};
        if (object->members.tail == -1)
            object->members.head = index;
        else
            data_.members.values()[static_cast<std::size_t>(object->members.tail)].next = index;
        object->members.tail = index;
        ++object->members.size;
        frame.pending_member = index;
        return true;
    }
    bool parse_error(std::size_t parser_position, const std::string&,
                     const Json::exception& error) override {
        event();
        // Never infer a source position from error.what(). Check the documented
        // one-based parser position, including the virtual byte after EOF.
        const auto n = static_cast<std::uint64_t>(data_.source.size());
        const auto position = static_cast<std::uint64_t>(parser_position);
        if (position != 0 && position - 1 <= n)
            context_.position = position - 1;
        const auto* syntax = dynamic_cast<const Json::parse_error*>(&error);
        const auto* range = dynamic_cast<const Json::out_of_range*>(&error);
        const std::string_view message{error.what()};
        // The pinned MSVC CRT silently drops its message on malloc failure.
        // Restrict this inference to parser-generated exception kinds, not an
        // arbitrary user/std exception with the same text.
        const auto kind = syntax != nullptr  ? detail::ParserFailureClass::Syntax
                          : range != nullptr ? detail::ParserFailureClass::NumberRange
                                             : detail::ParserFailureClass::Other;
        const auto classified = detail::classify_parser_failure(kind, error.id, message);
        context_.reject(classified.code, classified.reason);
    }

private:
    [[noreturn]] void impossible() {
        context_.reject(core::FailureCode::InvalidInput,
                        InputDiagnosticReason::UnexpectedException);
    }
    std::uint64_t current_event() {
        // Every callback publishes its event before storing a node or member.
        // Keep the checked observation local across other mutable helper calls.
        const auto ordinal = context_.event;
        if (!ordinal)
            impossible();
        return *ordinal;
    }
    void event() {
        context_.position.reset();
        context_.event = next_event_;
        next_event_ = add(next_event_, 1, context_);
    }
    detail::InputNode& at(core::Index value) noexcept {
        return data_.nodes.values()[static_cast<std::size_t>(value)];
    }
    void ensure_node() {
        const auto required = add(data_.node_count, 1, context_);
        if (required > limits_.maximum_nodes)
            context_.quota(InputQuota::Nodes, limits_.maximum_nodes, required);
        if (data_.node_count >= data_.nodes.size())
            impossible();
    }
    core::Index node(detail::NodeValue value) {
        ensure_node();
        const auto index = static_cast<core::Index>(data_.node_count++);
        at(index).value = value;
        at(index).event = current_event();
        if (depth_ == 0) {
            if (index != 0)
                impossible();
            return index;
        }
        auto& frame = frames_.values()[depth_ - 1];
        if (auto* array = std::get_if<detail::ArrayValue>(&at(frame.node).value)) {
            if (data_.edge_count >= data_.edges.size())
                impossible();
            const auto edge = static_cast<core::Index>(data_.edge_count++);
            data_.edges.values()[static_cast<std::size_t>(edge)].child = index;
            if (array->children.tail == -1)
                array->children.head = edge;
            else
                data_.edges.values()[static_cast<std::size_t>(array->children.tail)].next = edge;
            array->children.tail = edge;
            ++array->children.size;
        } else {
            if (!std::holds_alternative<detail::ObjectValue>(at(frame.node).value) ||
                frame.pending_member == -1)
                impossible();
            data_.members.values()[static_cast<std::size_t>(frame.pending_member)].child = index;
            frame.pending_member = -1;
        }
        return index;
    }
    template <class Container>
    void start(Container value) {
        const auto required = add(depth_, 1, context_);
        if (required > limits_.maximum_depth)
            context_.quota(InputQuota::Depth, limits_.maximum_depth, required);
        if (depth_ >= frames_.size())
            impossible();
        const auto index = node(value);
        frames_.values()[depth_++] = {index, -1};
    }
    void end(bool array) {
        if (depth_ == 0)
            impossible();
        const auto& frame = frames_.values()[depth_ - 1];
        if (frame.pending_member != -1 ||
            std::holds_alternative<detail::ArrayValue>(at(frame.node).value) != array)
            impossible();
        --depth_;
    }
    void check_text(std::string_view value, bool decimal) {
        const auto quota = decimal ? InputQuota::DecimalTokenBytes : InputQuota::DecodedStringBytes;
        const auto cap =
            decimal ? limits_.maximum_decimal_token_bytes : limits_.maximum_decoded_string_bytes;
        if (value.size() > cap)
            context_.quota(quota, cap, value.size());
        const auto required = add(data_.text_count, value.size(), context_);
        if (required > limits_.maximum_total_text_bytes)
            context_.quota(InputQuota::TotalTextBytes, limits_.maximum_total_text_bytes, required);
        if (required > data_.text.size())
            impossible();
    }
    detail::TextRange store(std::string_view value, bool decimal) {
        check_text(value, decimal);
        const detail::TextRange result{data_.text_count, value.size()};
        for (const char byte : value)
            data_.text.values()[data_.text_count++] = byte;
        return result;
    }
    static std::uint64_t key_hash(core::Index parent, std::string_view key) noexcept {
        // Hash is only an index accelerator: full parent/key bytes decide equality.
        std::uint64_t hash = UINT64_C(14695981039346656037) ^ static_cast<std::uint64_t>(parent);
        for (const char byte : key) {
            hash ^= static_cast<unsigned char>(byte);
            hash *= UINT64_C(1099511628211);
        }
        return hash;
    }
    InputDocumentData& data_;
    core::PlannedBuffer<detail::BuildFrame>& frames_;
    const InputLimits& limits_;
    Context& context_;
    std::size_t depth_ = 0;
    std::uint64_t next_event_ = 0;
};

void root_envelope(InputDocumentData& data, Context& context) {
    context.event = data.nodes.values()[0].event;
    context.position.reset();
    const auto* root = std::get_if<detail::ObjectValue>(&data.nodes.values()[0].value);
    if (root == nullptr)
        context.reject(core::FailureCode::InvalidInput, InputDiagnosticReason::RootType);
    constexpr std::array<std::string_view, 12> keys{
        "schema",   "units",     "model",  "flow",    "analysis", "observations",
        "numerics", "resources", "output", "initial", "restart",  "provenance"};
    const auto root_view = InputDocumentFactory::object(&data, 0);
    const auto schema = root_view.find("schema");
    if (schema) {
        context.event =
            data.nodes.values()[static_cast<std::size_t>(schema->ordinal().value())].event;
        const auto name = schema->string();
        if (!name)
            context.reject(core::FailureCode::InvalidInput, InputDiagnosticReason::SchemaType);
        if (*name != "uvlas.case/1")
            context.reject(core::FailureCode::Unsupported,
                           InputDiagnosticReason::UnsupportedSchema);
    }
    // Unknown keys are selected by decoded UTF-8 byte order, independently of
    // member order. Missing known keys use the fixed contract table order.
    core::Index unknown = -1;
    auto member = root->members.head;
    while (member != -1) {
        const auto key = member_key(data, member);
        if (std::find(keys.begin(), keys.end(), key) == keys.end() &&
            (unknown == -1 || key < member_key(data, unknown)))
            unknown = member;
        member = data.members.values()[static_cast<std::size_t>(member)].next;
    }
    if (unknown != -1) {
        context.event = data.members.values()[static_cast<std::size_t>(unknown)].event;
        context.reject(core::FailureCode::InvalidInput, InputDiagnosticReason::UnknownRootKey);
    }
    for (std::size_t key = 0; key < 9; ++key)
        if (!root_view.find(keys[key])) {
            context.event = data.nodes.values()[0].event;
            context.reject(core::FailureCode::InvalidInput, InputDiagnosticReason::MissingRootKey);
        }
}
} // namespace

core::Result<InputDocument> read_input_document(std::span<const std::byte> bytes, SourceId source,
                                                const InputLimits& limits,
                                                core::ResourcePlan& resources,
                                                std::optional<InputDiagnostic>& diagnostic,
                                                core::FailureDetails details) noexcept {
    diagnostic.reset();
    Context context{source,       diagnostic,   details,     std::nullopt,
                    std::nullopt, std::nullopt, std::nullopt};
    auto result = core::guard_entry(
        operation, core::FailureCode::InvalidInput, details, [&]() -> core::Result<InputDocument> {
            try {
                const auto n = static_cast<std::uint64_t>(bytes.size());
                if (n > limits.maximum_source_bytes)
                    context.quota(InputQuota::SourceBytes, limits.maximum_source_bytes, n);
                for (const auto cap :
                     {limits.maximum_source_bytes, limits.maximum_depth, limits.maximum_nodes,
                      limits.maximum_members, limits.maximum_decoded_string_bytes,
                      limits.maximum_decimal_token_bytes, limits.maximum_total_text_bytes})
                    (void)as_index(cap, context);
                if (limits.external_memory.mode != ParserAllocationMode::EstimateAdmission &&
                    limits.external_memory.mode != ParserAllocationMode::RequireProvenBound)
                    context.reject(core::FailureCode::InvalidInput,
                                   InputDiagnosticReason::InvalidLimits);
                context.plan = memory_plan(n, limits.maximum_depth, context);
                if (limits.external_memory.mode == ParserAllocationMode::RequireProvenBound)
                    context.reject(core::FailureCode::Unsupported,
                                   InputDiagnosticReason::UnprovenMemoryBound);
                if (context.plan->estimate_bytes() > limits.external_memory.maximum_bytes)
                    context.reject(core::FailureCode::ResourceExceeded,
                                   InputDiagnosticReason::ExternalMemoryLimit);
                const auto nodes = std::min(n, limits.maximum_nodes),
                           members = std::min(n, limits.maximum_members),
                           text_capacity = std::min(n, limits.maximum_total_text_bytes),
                           depth = std::min(n, limits.maximum_depth);
                const auto key_capacity = add(multiply(2, members, context), 1, context);
                (void)as_index(key_capacity, context);
                auto payload = static_cast<std::uint64_t>(sizeof(InputDocumentData));
                const auto account = [&](std::uint64_t count, std::size_t size) {
                    payload = add(payload, multiply(count, size, context), context);
                };
                account(n, sizeof(char));
                account(nodes, sizeof(detail::InputNode));
                account(nodes, sizeof(detail::ArrayEdge));
                account(members, sizeof(detail::ObjectMember));
                account(text_capacity, sizeof(char));
                account(key_capacity, sizeof(detail::KeySlot));
                account(depth, sizeof(detail::BuildFrame));
                if (payload > resources.remaining_bytes()) {
                    context.describe(InputDiagnosticReason::Allocation);
                    auto failure =
                        core::Failure{core::FailureCode::ResourceExceeded, operation, details};
                    const auto used = resources.byte_limit() - resources.remaining_bytes();
                    if (payload <= std::numeric_limits<std::uint64_t>::max() - used)
                        failure.details.resources =
                            core::ResourceUsage{core::ResourceKind::Bytes, used + payload, used};
                    return std::unexpected(failure);
                }
                auto owner = take(resources.allocate<InputDocumentData>(1, operation), context);
                auto& data = owner.values()[0];
                data.policy = limits.external_memory;
                data.parser_plan = context.plan;
                data.source =
                    take(resources.allocate<char>(as_index(n, context), operation), context);
                for (std::size_t i = 0; i < bytes.size(); ++i)
                    data.source.values()[i] = std::bit_cast<char>(bytes[i]);
                context.digest = detail::source_digest(std::as_bytes(data.source.values()));
                data.identity = SourceIdentity{source, n, *context.digest};
                if (bytes.empty())
                    context.reject(core::FailureCode::InvalidInput,
                                   InputDiagnosticReason::EmptyInput);
                const auto source_bytes = std::as_bytes(data.source.values());
                if (n >= 3 && source_bytes[0] == std::byte{0xef} &&
                    source_bytes[1] == std::byte{0xbb} && source_bytes[2] == std::byte{0xbf}) {
                    context.position = 0;
                    context.reject(core::FailureCode::InvalidInput,
                                   InputDiagnosticReason::ByteOrderMark);
                }
                for (std::size_t i = 0; i < source_bytes.size(); ++i)
                    if (source_bytes[i] == std::byte{0}) {
                        context.position = i;
                        context.reject(core::FailureCode::InvalidInput,
                                       InputDiagnosticReason::RawNull);
                    }
                data.nodes =
                    take(resources.allocate<detail::InputNode>(as_index(nodes, context), operation),
                         context);
                data.edges =
                    take(resources.allocate<detail::ArrayEdge>(as_index(nodes, context), operation),
                         context);
                data.members = take(
                    resources.allocate<detail::ObjectMember>(as_index(members, context), operation),
                    context);
                data.text = take(
                    resources.allocate<char>(as_index(text_capacity, context), operation), context);
                data.keys = take(
                    resources.allocate<detail::KeySlot>(as_index(key_capacity, context), operation),
                    context);
                auto frames = take(
                    resources.allocate<detail::BuildFrame>(as_index(depth, context), operation),
                    context);
                DocumentSax sax{data, frames, limits, context};
                auto adapter = nlohmann::detail::input_adapter(data.source.values().begin(),
                                                               data.source.values().end());
                nlohmann::detail::parser<Json, decltype(adapter)> parser(decltype(adapter){adapter},
                                                                         nullptr, true, false);
                if (!parser.sax_parse(&sax, true))
                    context.reject(core::FailureCode::InvalidInput,
                                   InputDiagnosticReason::JsonSyntax);
                root_envelope(data, context);
                diagnostic.reset();
                return InputDocumentFactory::document(std::move(owner));
            } catch (const Rejection& rejection) {
                return std::unexpected(rejection.failure);
            }
        });
    if (!result && !diagnostic) {
        context.describe(result.error().code == core::FailureCode::ResourceExceeded
                             ? InputDiagnosticReason::Allocation
                             : InputDiagnosticReason::UnexpectedException);
    }
    return result;
}
} // namespace uvlas::io
