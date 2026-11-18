#include "input_document_detail.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <memory>
#include <new>
#include <optional>
#include <string_view>
#if defined(_WIN32)
#include <malloc.h>
#endif

namespace {
// Process-local injection belongs only to this dedicated test executable.
constinit std::atomic<std::int64_t> remaining{-1};
constinit std::atomic<std::uint64_t> attempts{0};
// Replacement new injects only its standard allocation failure. Throwing
// arbitrary exceptions here can violate standard-library noexcept contracts;
// core/failure_test.cpp tests guard_entry's arbitrary-exception guarantee.
void allocation_attempt() {
    const auto next = remaining.load(std::memory_order_relaxed);
    if (next != -1) {
        attempts.fetch_add(1, std::memory_order_relaxed);
        if (next == 0) {
            // Disable before exception construction and stack unwinding.
            remaining.store(-1, std::memory_order_relaxed);
            throw std::bad_alloc{};
        }
        if (next > 0)
            remaining.fetch_sub(1, std::memory_order_relaxed);
    }
}
void* allocate(std::size_t count) {
    allocation_attempt();
    // NOLINTNEXTLINE(misc-const-correctness): operator new must return writable pointee storage.
    void* const value = std::malloc(count == 0 ? 1 : count);
    if (value == nullptr)
        throw std::bad_alloc{};
    return value;
}
void* allocate_aligned(std::size_t count, std::size_t alignment) {
    allocation_attempt();
    count = count == 0 ? 1 : count;
#if defined(_WIN32)
    // NOLINTNEXTLINE(misc-const-correctness): operator new must return writable storage.
    void* const value = _aligned_malloc(count, alignment);
#else
    if (count > std::numeric_limits<std::size_t>::max() - (alignment - 1))
        throw std::bad_alloc{};
    const auto rounded = ((count + alignment - 1) / alignment) * alignment;
    // NOLINTNEXTLINE(misc-const-correctness): operator new must return writable storage.
    void* const value = std::aligned_alloc(alignment, rounded);
#endif
    if (value == nullptr)
        throw std::bad_alloc{};
    return value;
}
void free_aligned(void* value) noexcept {
#if defined(_WIN32)
    _aligned_free(value);
#else
    std::free(value);
#endif
}
struct FaultScope {
    explicit FaultScope(std::int64_t count) noexcept {
        attempts.store(0, std::memory_order_relaxed);
        remaining.store(count, std::memory_order_relaxed);
    }
    ~FaultScope() {
        remaining.store(-1, std::memory_order_relaxed);
    }
    FaultScope(const FaultScope&) = delete;
    FaultScope& operator=(const FaultScope&) = delete;
};
} // namespace
void* operator new(std::size_t count) {
    return allocate(count);
}
void* operator new[](std::size_t count) {
    return allocate(count);
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
// Cover the forms used by library temporary buffers as well as user new
// expressions. Otherwise ASan's default nothrow new can pair with our free.
void* operator new(std::size_t count, const std::nothrow_t&) noexcept {
    try {
        return allocate(count);
    } catch (const std::bad_alloc&) {
        return nullptr;
    }
}
void* operator new[](std::size_t count, const std::nothrow_t&) noexcept {
    try {
        return allocate(count);
    } catch (const std::bad_alloc&) {
        return nullptr;
    }
}
void operator delete(void* value, const std::nothrow_t&) noexcept {
    std::free(value);
}
void operator delete[](void* value, const std::nothrow_t&) noexcept {
    std::free(value);
}
void* operator new(std::size_t count, std::align_val_t alignment) {
    return allocate_aligned(count, static_cast<std::size_t>(alignment));
}
void* operator new[](std::size_t count, std::align_val_t alignment) {
    return allocate_aligned(count, static_cast<std::size_t>(alignment));
}
void operator delete(void* value, std::align_val_t) noexcept {
    free_aligned(value);
}
void operator delete[](void* value, std::align_val_t) noexcept {
    free_aligned(value);
}
void operator delete(void* value, std::size_t, std::align_val_t) noexcept {
    free_aligned(value);
}
void operator delete[](void* value, std::size_t, std::align_val_t) noexcept {
    free_aligned(value);
}
void* operator new(std::size_t count, std::align_val_t alignment, const std::nothrow_t&) noexcept {
    try {
        return allocate_aligned(count, static_cast<std::size_t>(alignment));
    } catch (const std::bad_alloc&) {
        return nullptr;
    }
}
void* operator new[](std::size_t count, std::align_val_t alignment,
                     const std::nothrow_t&) noexcept {
    try {
        return allocate_aligned(count, static_cast<std::size_t>(alignment));
    } catch (const std::bad_alloc&) {
        return nullptr;
    }
}
void operator delete(void* value, std::align_val_t, const std::nothrow_t&) noexcept {
    free_aligned(value);
}
void operator delete[](void* value, std::align_val_t, const std::nothrow_t&) noexcept {
    free_aligned(value);
}

#include "case_schema_corpus.hpp"
#include "case_schema_snapshot.hpp"

#include <future>
#include <string>
#include <uvlas/io/case_schema.hpp>

namespace {
using namespace uvlas;
using namespace uvlas::io;
constexpr std::uint64_t schema_budget = 256ULL * 1024 * 1024;
constexpr std::size_t quota_count = static_cast<std::size_t>(SchemaQuota::Count);
using Quotas = std::array<std::uint64_t, quota_count>;
Quotas generous_quotas() noexcept {
    Quotas values{};
    values.fill(1000000);
    values[static_cast<std::size_t>(SchemaQuota::TypedBytes)] = schema_budget;
    values[static_cast<std::size_t>(SchemaQuota::ScratchBytes)] = schema_budget;
    values[static_cast<std::size_t>(SchemaQuota::Work)] = 1000000000;
    values[static_cast<std::size_t>(SchemaQuota::ExactDigits)] = 4096;
    return values;
}
core::Result<InputDocument> read_schema_source(std::string_view source, core::ResourcePlan& plan) {
    const InputLimits syntax{
        1000000,
        256,
        100000,
        100000,
        1000000,
        1000000,
        1000000,
        {ParserAllocationMode::EstimateAdmission, std::numeric_limits<std::uint64_t>::max()}};
    std::optional<InputDiagnostic> diagnostic;
    return read_input_document(std::as_bytes(std::span{source.data(), source.size()}), SourceId{7},
                               syntax, plan, diagnostic);
}
struct SchemaSnapshot {
    SourceIdentity identity;
    std::uint64_t geometry_count;
    std::uint64_t shell_count;
    std::uint64_t beam_count;
    std::uint64_t observation_count;
    std::uint64_t owned_bytes;
    std::uint64_t scratch_bytes;
    std::uint64_t cumulative_bytes;
    std::uint64_t work;
    uvlas::test::schema_snapshot::Bytes payload;
    bool operator==(const SchemaSnapshot&) const = default;
};
core::Result<SchemaSnapshot> snapshot(const SchemaCase& value) {
    const auto root = value.root();
    const auto identity = value.source_identity();
    const auto memory = value.memory_report();
    if (!root)
        return std::unexpected(root.error());
    if (!identity)
        return std::unexpected(identity.error());
    if (!memory)
        return std::unexpected(memory.error());
    const auto model = root->model.record();
    return SchemaSnapshot{*identity,
                          model.geometry.size(),
                          model.shells.size(),
                          model.beams.size(),
                          root->observations.size(),
                          memory->owned_bytes,
                          memory->scratch_bytes,
                          memory->cumulative_bytes,
                          memory->work,
                          uvlas::test::schema_snapshot::serialize(*root)};
}
bool same_snapshot(const SchemaCase& value, const SchemaSnapshot& expected) {
    const auto observed = snapshot(value);
    return observed && *observed == expected;
}
} // namespace

TEST_CASE("Schema allocation failures preserve source old result and successful charges",
          "[negative][property][IO4-01][N4-02][N4-08]") {
    for (const auto& example : uvlas::test::case_schema_corpus::roots) {
        INFO(example.name);
        core::ResourcePlan input_plan{schema_budget};
        const auto document = read_schema_source(example.json, input_plan);
        REQUIRE(document);
        const auto view = document->view();
        REQUIRE(view);
        const auto source_identity = view->identity;
        const std::string original{example.json};
        core::ResourcePlan old_plan{schema_budget};
        SchemaDiagnostic diagnostic;
        const auto old =
            validate_case_schema(*document, SchemaLimits{generous_quotas()}, old_plan, &diagnostic);
        REQUIRE(old);
        const auto old_snapshot = snapshot(*old);
        REQUIRE(old_snapshot);
        core::ResourcePlan control_plan{schema_budget};
        const auto control = [&] {
            const FaultScope count{-2};
            return validate_case_schema(*document, SchemaLimits{generous_quotas()}, control_plan,
                                        &diagnostic);
        }();
        REQUIRE(control);
        const auto sites = attempts.load(std::memory_order_relaxed);
        REQUIRE(sites > 10);
        REQUIRE(sites < 10000);
        bool retained_charge = false;
        for (std::uint64_t site = 0; site < sites; ++site) {
            INFO(site);
            core::ResourcePlan plan{schema_budget};
            const auto prior = plan.allocate<std::byte>(3, core::Operation::PrepareModel);
            REQUIRE(prior);
            const auto failed = [&] {
                const FaultScope fault{static_cast<std::int64_t>(site)};
                return validate_case_schema(*document, SchemaLimits{generous_quotas()}, plan,
                                            &diagnostic);
            }();
            REQUIRE_FALSE(failed);
            CHECK(failed.error().code == core::FailureCode::ResourceExceeded);
            CHECK(failed.error().operation == core::Operation::PrepareModel);
            const auto detail = diagnostic.detail();
            if (!detail) {
                FAIL("Expected owned schema diagnostic is absent");
                return;
            }
            CHECK(detail->reason == SchemaReason::Allocation);
            CHECK(detail->source == source_identity);
            CHECK(plan.allocated_bytes() >= 3);
            retained_charge = retained_charge || plan.allocated_bytes() > 3;
            CHECK(same_snapshot(*old, *old_snapshot));
            const auto unchanged = document->view();
            REQUIRE(unchanged);
            CHECK(unchanged->identity == source_identity);
            CHECK(std::equal(unchanged->source_bytes.begin(), unchanged->source_bytes.end(),
                             std::as_bytes(std::span{original.data(), original.size()}).begin()));
        }
        CHECK(retained_charge);
        core::ResourcePlan recovery_plan{schema_budget};
        CHECK(validate_case_schema(*document, SchemaLimits{generous_quotas()}, recovery_plan,
                                   &diagnostic));
        CHECK_FALSE(diagnostic.detail());
        std::printf("schema root=%.*s actual new fault sites=%llu\n",
                    static_cast<int>(example.name.size()), example.name.data(),
                    static_cast<unsigned long long>(sites));
    }
}

TEST_CASE("Schema source and old result survive independent simultaneous calls",
          "[property][IO4-01][N4-08]") {
    for (const auto& example : uvlas::test::case_schema_corpus::roots) {
        INFO(example.name);
        core::ResourcePlan input_plan{schema_budget};
        const auto document = read_schema_source(example.json, input_plan);
        REQUIRE(document);
        const auto worker = [&document]() -> core::Result<SchemaCase> {
            core::ResourcePlan plan{schema_budget};
            SchemaDiagnostic diagnostic;
            return validate_case_schema(*document, SchemaLimits{generous_quotas()}, plan,
                                        &diagnostic);
        };
        auto first_job = std::async(std::launch::async, worker);
        auto second_job = std::async(std::launch::async, worker);
        auto first = first_job.get();
        auto second = second_job.get();
        REQUIRE(first);
        REQUIRE(second);
        const auto first_snapshot = snapshot(*first);
        const auto second_snapshot = snapshot(*second);
        REQUIRE(first_snapshot);
        REQUIRE(second_snapshot);
        CHECK(*first_snapshot == *second_snapshot);
        core::ResourcePlan no_slot_plan{schema_budget};
        const auto no_slot =
            validate_case_schema(*document, SchemaLimits{generous_quotas()}, no_slot_plan, nullptr);
        REQUIRE(no_slot);
        CHECK(same_snapshot(*no_slot, *first_snapshot));
        auto moved = std::move(*first);
        CHECK_FALSE(first->root());
        CHECK(same_snapshot(moved, *second_snapshot));
    }
}

TEST_CASE("Every schema quota has a finite monotone acceptance boundary",
          "[property][negative][IO4-01][N4-02]") {
    for (const auto& example : uvlas::test::case_schema_corpus::roots) {
        INFO(example.name);
        core::ResourcePlan input_plan{schema_budget};
        const auto document = read_schema_source(example.json, input_plan);
        REQUIRE(document);
        const auto quotas = generous_quotas();
        for (std::size_t index = 0; index < quota_count; ++index) {
            INFO(index);
            const auto accepts = [&](std::uint64_t limit) {
                auto requested = quotas;
                requested[index] = limit;
                core::ResourcePlan plan{schema_budget};
                SchemaDiagnostic diagnostic;
                const auto result =
                    validate_case_schema(*document, SchemaLimits{requested}, plan, &diagnostic);
                if (!result) {
                    CHECK(result.error().code == core::FailureCode::ResourceExceeded);
                    const auto detail = diagnostic.detail();
                    if (!detail) {
                        FAIL("Expected owned schema diagnostic is absent");
                        return false;
                    }
                    CHECK((detail->reason == SchemaReason::QuotaExceeded ||
                           detail->reason == SchemaReason::Overflow));
                    if (detail->reason == SchemaReason::QuotaExceeded)
                        CHECK(detail->quota == static_cast<SchemaQuota>(index));
                }
                return result.has_value();
            };
            REQUIRE(accepts(quotas[index]));
            std::uint64_t lower = 0;
            auto upper = quotas[index];
            while (lower < upper) {
                const auto middle = lower + (upper - lower) / 2;
                if (accepts(middle))
                    upper = middle;
                else
                    lower = middle + 1;
            }
            CHECK(accepts(lower));
            if (lower > 0)
                CHECK_FALSE(accepts(lower - 1));
            // An unused quota may legitimately have a minimum of zero.
            std::printf("schema root=%.*s quota=%zu minimum=%llu\n",
                        static_cast<int>(example.name.size()), example.name.data(), index,
                        static_cast<unsigned long long>(lower));
        }
    }
}

#include <vector>

namespace {
struct SourceDemands {
    std::array<std::uint64_t, 7> original{};
    std::uint64_t path_bytes = 0;
};
SourceDemands source_demands(const InputDocument& document) {
    const auto* source = uvlas::io::detail::InputDocumentFactory::storage(document);
    SourceDemands result;
    if (source == nullptr)
        return result;
    result.original[0] = source->source.size();
    result.original[2] = source->node_count;
    result.original[3] = source->member_count;
    result.original[6] = source->text_count;
    std::vector<std::uint64_t> depths(source->node_count, 0);
    std::vector<std::string> paths(source->node_count);
    depths[0] = 1;
    const auto escape = [](std::string_view key) {
        std::string escaped;
        for (const char ch : key) {
            if (ch == '~')
                escaped += "~0";
            else if (ch == '/')
                escaped += "~1";
            else
                escaped += ch;
        }
        return escaped;
    };
    for (std::size_t index = 0; index < source->node_count; ++index) {
        result.path_bytes =
            std::max(result.path_bytes, static_cast<std::uint64_t>(paths[index].size()));
        const auto& value = source->nodes.values()[index].value;
        if (const auto* text = std::get_if<uvlas::io::detail::StringValue>(&value))
            result.original[4] =
                std::max(result.original[4], static_cast<std::uint64_t>(text->text.size));
        if (const auto* decimal = std::get_if<uvlas::io::detail::DecimalValue>(&value))
            result.original[5] =
                std::max(result.original[5], static_cast<std::uint64_t>(decimal->text.size));
        if (const auto* object = std::get_if<uvlas::io::detail::ObjectValue>(&value)) {
            result.original[1] = std::max(result.original[1], depths[index]);
            auto member = object->members.head;
            while (member >= 0) {
                const auto& entry = source->members.values()[static_cast<std::size_t>(member)];
                const auto child = static_cast<std::size_t>(entry.child);
                const std::string_view key{source->text.values().data() + entry.key.offset,
                                           entry.key.size};
                result.original[4] =
                    std::max(result.original[4], static_cast<std::uint64_t>(key.size()));
                depths[child] = depths[index] + 1;
                paths[child] = paths[index] + "/" + escape(key);
                member = entry.next;
            }
        }
        if (const auto* array = std::get_if<uvlas::io::detail::ArrayValue>(&value)) {
            result.original[1] = std::max(result.original[1], depths[index]);
            auto edge = array->children.head;
            std::uint64_t ordinal = 0;
            while (edge >= 0) {
                const auto& entry = source->edges.values()[static_cast<std::size_t>(edge)];
                const auto child = static_cast<std::size_t>(entry.child);
                depths[child] = depths[index] + 1;
                paths[child] = paths[index] + "/" + std::to_string(ordinal++);
                edge = entry.next;
            }
        }
    }
    return result;
}
} // namespace

TEST_CASE("Schema admission matches independent source counts and explicit payload reports",
          "[verification][negative][IO4-01][N4-02]") {
    for (const auto& example : uvlas::test::case_schema_corpus::roots) {
        INFO(example.name);
        core::ResourcePlan input_plan{schema_budget};
        const auto document = read_schema_source(example.json, input_plan);
        REQUIRE(document);
        const auto original = source_demands(*document);
        core::ResourcePlan control_plan{schema_budget};
        const auto control =
            validate_case_schema(*document, SchemaLimits{generous_quotas()}, control_plan, nullptr);
        REQUIRE(control);
        const auto report = control->memory_report();
        REQUIRE(report);
        CHECK(report->owned_bytes <= report->planned_owned_bytes);
        CHECK(report->scratch_bytes <= report->planned_scratch_bytes);
        CHECK(report->cumulative_bytes == control_plan.allocated_bytes());
        CHECK(report->cumulative_bytes <= report->admission_bytes);
        CHECK(report->admission_bytes ==
              report->planned_owned_bytes + report->planned_scratch_bytes + original.path_bytes);
        std::array<std::uint64_t, quota_count> demands{};
        std::copy(original.original.begin(), original.original.end(), demands.begin());
        demands[static_cast<std::size_t>(SchemaQuota::TypedBytes)] = report->planned_owned_bytes;
        demands[static_cast<std::size_t>(SchemaQuota::ScratchBytes)] =
            report->planned_scratch_bytes;
        demands[static_cast<std::size_t>(SchemaQuota::Work)] = report->work;
        demands[static_cast<std::size_t>(SchemaQuota::DiagnosticPathBytes)] = original.path_bytes;
        for (std::size_t index = 0; index < quota_count; ++index) {
            if (demands[index] == 0)
                continue;
            INFO(index);
            auto exact = generous_quotas();
            exact[index] = demands[index];
            core::ResourcePlan exact_plan{schema_budget};
            SchemaDiagnostic diagnostic;
            CHECK(validate_case_schema(*document, SchemaLimits{exact}, exact_plan, &diagnostic));
            CHECK_FALSE(diagnostic.detail());
            --exact[index];
            core::ResourcePlan shortage_plan{schema_budget};
            const auto failed =
                validate_case_schema(*document, SchemaLimits{exact}, shortage_plan, &diagnostic);
            REQUIRE_FALSE(failed);
            CHECK(failed.error().code == core::FailureCode::ResourceExceeded);
            const auto detail = diagnostic.detail();
            if (!detail) {
                FAIL("Expected owned schema diagnostic is absent");
                return;
            }
            CHECK(detail->quota == static_cast<SchemaQuota>(index));
        }
        // Admission includes all successful allocations plus the potential failure path.
        // A previously charged allocation cannot reuse that reserved capacity.
        constexpr std::uint64_t prior_bytes = 11;
        for (const auto shortage : {std::uint64_t{0}, std::uint64_t{1}}) {
            core::ResourcePlan plan{report->admission_bytes + prior_bytes - shortage};
            const auto prior = plan.allocate<std::byte>(prior_bytes, core::Operation::PrepareModel);
            REQUIRE(prior);
            SchemaDiagnostic diagnostic;
            const auto result =
                validate_case_schema(*document, SchemaLimits{generous_quotas()}, plan, &diagnostic);
            CHECK(result.has_value() == (shortage == 0));
            CHECK(plan.allocated_bytes() >= prior_bytes);
            if (!result)
                CHECK(result.error().code == core::FailureCode::ResourceExceeded);
        }
    }
}

TEST_CASE("Schema path failures retain the original schema rejection and own their context",
          "[property][negative][IO4-01][N4-08]") {
    constexpr std::string_view invalid =
        R"({"schema":"uvlas.case/1","units":null,"model":null,"flow":null,"analysis":null,"observations":null,"numerics":null,"resources":null,"output":null})";
    core::ResourcePlan input_plan{schema_budget};
    const auto document = read_schema_source(invalid, input_plan);
    REQUIRE(document);
    SchemaDiagnostic diagnostic;
    core::ResourcePlan control_plan{schema_budget};
    const auto control = [&] {
        const FaultScope count{-2};
        return validate_case_schema(*document, SchemaLimits{generous_quotas()}, control_plan,
                                    &diagnostic);
    }();
    REQUIRE_FALSE(control);
    const auto sites = attempts.load(std::memory_order_relaxed);
    const auto original = diagnostic.detail();
    if (!original) {
        FAIL("Expected owned schema diagnostic is absent");
        return;
    }
    CHECK(original->reason == SchemaReason::WrongKind);
    const auto original_path = diagnostic.path();
    if (!original_path) {
        FAIL("Expected owned schema diagnostic is absent");
        return;
    }
    CHECK(*original_path == "/units");
    REQUIRE(sites >= 2);
    core::ResourcePlan failed_path_plan{schema_budget};
    const auto failed_path = [&] {
        const FaultScope fault{static_cast<std::int64_t>(sites - 1)};
        return validate_case_schema(*document, SchemaLimits{generous_quotas()}, failed_path_plan,
                                    &diagnostic);
    }();
    REQUIRE_FALSE(failed_path);
    CHECK(failed_path.error().code == control.error().code);
    const auto path_detail = diagnostic.detail();
    if (!path_detail) {
        FAIL("Expected owned schema diagnostic is absent");
        return;
    }
    CHECK(path_detail->reason == original->reason);
    CHECK(path_detail->source == original->source);
    CHECK(path_detail->node == original->node);
    CHECK(path_detail->path_state == SchemaPathState::Allocation);
    CHECK_FALSE(diagnostic.path());
    auto no_path_limits = generous_quotas();
    no_path_limits[static_cast<std::size_t>(SchemaQuota::DiagnosticPathBytes)] = 0;
    for (const bool with_slot : {false, true}) {
        core::ResourcePlan plan{schema_budget};
        const auto result = validate_case_schema(*document, SchemaLimits{no_path_limits}, plan,
                                                 with_slot ? &diagnostic : nullptr);
        REQUIRE_FALSE(result);
        CHECK(result.error().code == control.error().code);
        if (with_slot) {
            const auto detail = diagnostic.detail();
            if (!detail) {
                FAIL("Expected owned schema diagnostic is absent");
                return;
            }
            CHECK(detail->reason == original->reason);
            CHECK(detail->path_state == SchemaPathState::PathQuota);
            CHECK_FALSE(diagnostic.path());
        }
    }
}

TEST_CASE("Every case input declaration can lower but cannot enlarge its host quota",
          "[property][negative][IO4-01][N4-02]") {
    constexpr std::array names{"maximum_source_bytes",
                               "maximum_depth",
                               "maximum_nodes",
                               "maximum_members",
                               "maximum_decoded_string_bytes",
                               "maximum_decimal_token_bytes",
                               "maximum_total_text_bytes",
                               "maximum_typed_bytes",
                               "maximum_entities",
                               "maximum_references",
                               "maximum_scalars",
                               "maximum_polynomial_terms",
                               "maximum_table_samples",
                               "maximum_expression_nodes",
                               "maximum_expression_edges",
                               "maximum_expression_depth",
                               "maximum_expression_power",
                               "maximum_expression_points",
                               "maximum_diagnostic_path_bytes"};
    for (const auto& example : uvlas::test::case_schema_corpus::roots) {
        INFO(example.name);
        for (std::size_t index = 0; index < names.size(); ++index) {
            INFO(names[index]);
            const std::string key = std::string{"\""} + names[index] + "\"";
            const auto key_position = example.json.find(key);
            REQUIRE(key_position != std::string_view::npos);
            REQUIRE(example.json.find(key, key_position + key.size()) == std::string_view::npos);
            auto begin = example.json.find(':', key_position + key.size());
            REQUIRE(begin != std::string_view::npos);
            ++begin;
            while (begin < example.json.size() &&
                   (example.json[begin] == ' ' || example.json[begin] == '\n' ||
                    example.json[begin] == '\r' || example.json[begin] == '\t'))
                ++begin;
            auto end = begin;
            while (end < example.json.size() && example.json[end] >= '0' &&
                   example.json[end] <= '9')
                ++end;
            REQUIRE(end > begin);
            const auto run = [&](std::string_view declaration, bool host_zero) {
                std::string source{example.json};
                source.replace(begin, end - begin, declaration);
                core::ResourcePlan input_plan{schema_budget};
                auto document = read_schema_source(source, input_plan);
                if (!document) {
                    FAIL_CHECK("Changing a native quota must preserve syntax");
                    return false;
                }
                auto requested = generous_quotas();
                if (host_zero)
                    requested[index] = 0;
                core::ResourcePlan plan{schema_budget};
                SchemaDiagnostic diagnostic;
                const auto result =
                    validate_case_schema(*document, SchemaLimits{requested}, plan, &diagnostic);
                if (!result) {
                    CHECK(result.error().code == core::FailureCode::ResourceExceeded);
                    const auto detail = diagnostic.detail();
                    if (!detail) {
                        FAIL_CHECK("Quota rejection must own its cause");
                        return false;
                    }
                    CHECK(detail->quota == static_cast<SchemaQuota>(index));
                }
                return result.has_value();
            };
            const auto case_zero = run("0", false);
            const auto host_zero = run("18446744073709551615", true);
            CHECK(case_zero == host_zero);
            // Source/shape/owned/path quotas are necessarily nonzero in every full case.
            if (index <= static_cast<std::size_t>(SchemaQuota::TypedBytes) ||
                index == static_cast<std::size_t>(SchemaQuota::DiagnosticPathBytes))
                CHECK_FALSE(case_zero);
        }
    }
}
