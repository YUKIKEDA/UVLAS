// Allocation replacement follows the existing typed-schema resource executable.

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

#include "case_schema_test_adapter.hpp"
#include "si_case_ownership_snapshot.hpp"

#include <future>
#include <utility>
#include <uvlas/io/si_case.hpp>
#include <vector>

namespace {
using namespace uvlas;
using namespace uvlas::io;
using test::schema_adapter::root;
using test::schema_adapter::string;
using test::si_ownership_snapshot::serialize;

constexpr std::uint64_t budget = std::uint64_t{1} << 50U;
using Quotas = std::array<std::uint64_t, static_cast<std::size_t>(SiQuota::Count)>;

constexpr std::size_t position(SiQuota quota) {
    return static_cast<std::size_t>(quota);
}

Quotas generous() {
    Quotas values{};
    values.fill(1'000'000);
    for (const auto quota :
         {SiQuota::OwnedBytes, SiQuota::ScratchBytes, SiQuota::CumulativeBytes, SiQuota::Work})
        values[position(quota)] = budget;
    return values;
}

SiLimits limits(const Quotas& values, std::stop_token stop = {}) {
    return {values,
            {100000, 100000, 100000, 64, 64, 1000000, 2000000,
             structure::ExactEstimateAdmission{budget}, structure::ExactTimingMode::Disabled},
            std::move(stop)};
}

template <class T>
T take(core::Result<T> result) {
    REQUIRE(result);
    return std::move(*result);
}

SchemaCase source(std::string_view reference_length = {}) {
    auto fixture = root();
    fixture.field("units").set("length", string("mm"));
    if (!reference_length.empty())
        fixture.field("model")
            .field("geometry")
            .elements.at(0)
            .set("reference_length", test::schema_adapter::parse(reference_length));
    const auto text = fixture.wire();
    const InputLimits input{
        1'000'000, 256,       100'000,   100'000,
        1'000'000, 1'000'000, 1'000'000, {ParserAllocationMode::EstimateAdmission, budget}};
    core::ResourcePlan plan{budget};
    std::optional<InputDiagnostic> diagnostic;
    const auto document =
        take(read_input_document(std::as_bytes(std::span{text.data(), text.size()}), SourceId{718},
                                 input, plan, diagnostic));
    std::array<std::uint64_t, static_cast<std::size_t>(SchemaQuota::Count)> caps{};
    caps.fill(1'000'000'000);
    return take(validate_case_schema(document, SchemaLimits{caps}, plan, nullptr));
}

std::vector<SiRequirementRow> requirements(const SiCatalog& catalog) {
    std::vector<SiRequirementRow> rows;
    for (std::uint64_t i = 0; i < take(catalog.size()); ++i) {
        const auto entry = take(catalog.entry(i));
        rows.push_back({entry.key,
                        {SiDecimalInput{"1"}, entry.recipe.dimensions},
                        SiDecimalInput{"1e-25"},
                        SiDecimalInput{"1e-25"}});
    }
    return rows;
}

template <class Call>
void check_accounting(const Call& call) {
    core::ResourcePlan generous_plan{budget};
    const auto result = call(limits(generous()), generous_plan, nullptr);
    REQUIRE(result);
    const auto report = take(result->memory_report());
    CHECK(report.cumulative_bytes == generous_plan.allocated_bytes());
    CHECK(report.owned_bytes <= report.planned_owned_bytes);
    CHECK(report.scratch_bytes <= report.planned_scratch_bytes);
    CHECK(report.cumulative_bytes <= report.admission_bytes);
    CHECK(report.work <= report.planned_work);
    for (const auto& [quota, required] :
         std::array{std::pair{SiQuota::OwnedBytes, report.planned_owned_bytes},
                    std::pair{SiQuota::ScratchBytes, report.planned_scratch_bytes},
                    std::pair{SiQuota::CumulativeBytes, report.admission_bytes},
                    std::pair{SiQuota::Work, report.planned_work}}) {
        INFO(position(quota));
        if (required == 0) {
            CHECK(quota == SiQuota::ScratchBytes);
            continue;
        }
        auto caps = generous();
        caps[position(quota)] = required;
        core::ResourcePlan equality{budget};
        CHECK(call(limits(caps), equality, nullptr));
        caps[position(quota)] = required - 1;
        core::ResourcePlan shortage{budget};
        SiDiagnostic detail;
        const auto rejected = call(limits(caps), shortage, &detail);
        REQUIRE_FALSE(rejected);
        CHECK(rejected.error().code == core::FailureCode::ResourceExceeded);
        const auto detail_view = detail.detail();
        REQUIRE((detail_view && detail_view->quota == quota));
    }
    core::ResourcePlan bytes_equal{report.admission_bytes};
    CHECK(call(limits(generous()), bytes_equal, nullptr));
    REQUIRE(report.admission_bytes > 0);
    core::ResourcePlan bytes_short{report.admission_bytes - 1};
    const auto failed = call(limits(generous()), bytes_short, nullptr);
    REQUIRE_FALSE(failed);
    CHECK(failed.error().code == core::FailureCode::ResourceExceeded);
    CHECK(bytes_short.allocated_bytes() == 0);
}

template <class Call>
void check_source_quotas(const Call& call, std::span<const SiQuota> quotas) {
    for (const auto quota : quotas) {
        INFO(position(quota));
        // The small published fixture bounds the range. Locate the exact accepted
        // count independently of the success report and test both sides explicitly.
        const auto attempt = [&](std::uint64_t maximum) {
            auto configured = generous();
            configured[position(quota)] = maximum;
            core::ResourcePlan plan{budget};
            return call(limits(configured), plan, nullptr);
        };
        REQUIRE(attempt(1'000'000));
        const auto zero = attempt(0);
        REQUIRE_FALSE(zero);
        REQUIRE(zero.error().code == core::FailureCode::ResourceExceeded);
        std::uint64_t low = 0, high = 1'000'000;
        while (high - low > 1) {
            const auto middle = low + (high - low) / 2;
            const auto result = attempt(middle);
            if (result)
                high = middle;
            else {
                REQUIRE(result.error().code == core::FailureCode::ResourceExceeded);
                low = middle;
            }
        }
        REQUIRE(attempt(high));
        REQUIRE_FALSE(attempt(high - 1));
    }
}

template <class Call, class CheckParent>
void fail_allocations(const Call& call, const CheckParent& check_parent) {
    const auto admitted = limits(generous());
    core::ResourcePlan plan{budget};
    SiDiagnostic diagnostic;
    const auto control = [&] {
        const FaultScope count{-2};
        return call(admitted, plan, &diagnostic);
    }();
    REQUIRE(control);
    const auto sites = attempts.load(std::memory_order_relaxed);
    REQUIRE(sites > 0);
    REQUIRE(sites < 10'000);
    bool paid_prefix = false;
    for (std::uint64_t site = 0; site < sites; ++site) {
        INFO(site);
        const auto before = plan.allocated_bytes();
        const auto failed = [&] {
            const FaultScope inject{static_cast<std::int64_t>(site)};
            return call(admitted, plan, &diagnostic);
        }();
        REQUIRE_FALSE(failed);
        CHECK(failed.error().code == core::FailureCode::ResourceExceeded);
        CHECK(failed.error().operation == core::Operation::PrepareModel);
        CHECK(plan.allocated_bytes() >= before);
        paid_prefix = paid_prefix || plan.allocated_bytes() > before;
        check_parent(diagnostic);
    }
    CHECK(paid_prefix);
    REQUIRE(call(admitted, plan, &diagnostic));
    CHECK_FALSE(diagnostic.detail());
}
} // namespace

TEST_CASE("SI owns prior results through every observed allocation failure",
          "[negative][verification][IO4-01][N4-02]") {
    const auto schema = source();
    core::ResourcePlan plan{budget};
    const auto admitted = limits(generous());
    const auto catalog = take(make_si_catalog(schema, admitted, plan, nullptr));
    const auto rows = requirements(catalog);
    const auto request = take(make_si_request(catalog, rows, admitted, plan, nullptr));
    const auto old = take(convert_case_to_si(schema, catalog, request, admitted, plan, nullptr));
    const auto unchanged = serialize(old);
    const auto compare = [&](const SiDiagnostic&) { CHECK(serialize(old) == unchanged); };
    fail_allocations(
        [&](const SiLimits& value, core::ResourcePlan& usage, SiDiagnostic* diagnostic) {
            return make_si_catalog(schema, value, usage, diagnostic);
        },
        compare);
    fail_allocations(
        [&](const SiLimits& value, core::ResourcePlan& usage, SiDiagnostic* diagnostic) {
            return make_si_request(catalog, rows, value, usage, diagnostic);
        },
        compare);
    fail_allocations(
        [&](const SiLimits& value, core::ResourcePlan& usage, SiDiagnostic* diagnostic) {
            return convert_case_to_si(schema, catalog, request, value, usage, diagnostic);
        },
        compare);
}

TEST_CASE("All SI quotas admit their boundary and reject one unit before it",
          "[negative][verification][IO4-01][N4-02]") {
    const auto schema = source();
    core::ResourcePlan plan{budget};
    const auto admitted = limits(generous());
    const auto catalog = take(make_si_catalog(schema, admitted, plan, nullptr));
    const auto rows = requirements(catalog);
    const auto request = take(make_si_request(catalog, rows, admitted, plan, nullptr));
    const auto catalog_call = [&](const SiLimits& value, core::ResourcePlan& usage,
                                  SiDiagnostic* diagnostic) {
        return make_si_catalog(schema, value, usage, diagnostic);
    };
    const auto request_call = [&](const SiLimits& value, core::ResourcePlan& usage,
                                  SiDiagnostic* diagnostic) {
        return make_si_request(catalog, rows, value, usage, diagnostic);
    };
    const auto convert_call = [&](const SiLimits& value, core::ResourcePlan& usage,
                                  SiDiagnostic* diagnostic) {
        return convert_case_to_si(schema, catalog, request, value, usage, diagnostic);
    };
    const std::array source_caps{
        SiQuota::SourceBytes,       SiQuota::InputNodes,   SiQuota::InputMembers,
        SiQuota::InputDepth,        SiQuota::SourceDigits, SiQuota::Quantities,
        SiQuota::ExponentMagnitude, SiQuota::UnitPower,    SiQuota::DiagnosticPathBytes};
    check_source_quotas(catalog_call, source_caps);
    const std::array request_caps{SiQuota::RequestTokenBytes, SiQuota::RequestDigits};
    check_source_quotas(request_call, request_caps);
    check_accounting(catalog_call);
    check_accounting(request_call);
    check_accounting(convert_call);
}

TEST_CASE("SI cancellation and zero diagnostic capacity preserve their original failures",
          "[negative][IO4-01][N4-02]") {
    const auto schema = source();
    core::ResourcePlan plan{budget};
    const auto admitted = limits(generous());
    const auto catalog = take(make_si_catalog(schema, admitted, plan, nullptr));
    auto rows = requirements(catalog);
    const auto request = take(make_si_request(catalog, rows, admitted, plan, nullptr));
    const auto old = take(convert_case_to_si(schema, catalog, request, admitted, plan, nullptr));
    const auto unchanged = serialize(old);
    std::stop_source stop;
    REQUIRE(stop.request_stop());
    const auto stopped = limits(generous(), stop.get_token());
    const auto before = plan.allocated_bytes();
    SiDiagnostic detail;
    const auto stopped_catalog = make_si_catalog(schema, stopped, plan, &detail);
    const auto stopped_request = make_si_request(catalog, rows, stopped, plan, &detail);
    const auto stopped_case = convert_case_to_si(schema, catalog, request, stopped, plan, &detail);
    REQUIRE_FALSE(stopped_catalog);
    REQUIRE_FALSE(stopped_request);
    REQUIRE_FALSE(stopped_case);
    for (const auto code :
         {stopped_catalog.error().code, stopped_request.error().code, stopped_case.error().code})
        CHECK(code == core::FailureCode::Cancelled);
    CHECK(plan.allocated_bytes() == before);

    rows.pop_back();
    auto caps = generous();
    caps[position(SiQuota::CumulativeBytes)] = 0;
    const auto failed = make_si_request(catalog, rows, limits(caps), plan, &detail);
    REQUIRE_FALSE(failed);
    CHECK(failed.error().code == core::FailureCode::InvalidInput);
    const auto missing_detail = detail.detail();
    REQUIRE((missing_detail && missing_detail->reason == SiReason::MissingRequirement));
    CHECK_FALSE(detail.path());
    CHECK(plan.allocated_bytes() == before);

    // SourceDigits rejects at a real parsed field before admission. A path is
    // available without the cumulative cap, so this exercises the diagnostic's
    // own check rather than an error which has no location to begin with.
    caps[position(SiQuota::SourceDigits)] = 0;
    const auto capped = make_si_catalog(schema, limits(caps), plan, &detail);
    REQUIRE_FALSE(capped);
    const auto capped_detail = detail.detail();
    REQUIRE((capped_detail && capped_detail->quota == SiQuota::SourceDigits));
    CHECK(capped.error().code == core::FailureCode::ResourceExceeded);
    CHECK_FALSE(detail.path());
    CHECK(plan.allocated_bytes() == before);
    caps[position(SiQuota::CumulativeBytes)] = budget;
    const auto located = make_si_catalog(schema, limits(caps), plan, &detail);
    REQUIRE_FALSE(located);
    CHECK(located.error().code == capped.error().code);
    const auto located_path = detail.path();
    REQUIRE((located_path && !located_path->empty()));
    CHECK(serialize(old) == unchanged);
}

TEST_CASE("Independent SI calls retain their exact payload under concurrent execution",
          "[property][verification][IO4-01][N4-08]") {
    const auto schema = source();
    const auto admitted = limits(generous());
    core::ResourcePlan plan{budget};
    const auto catalog = take(make_si_catalog(schema, admitted, plan, nullptr));
    const auto rows = requirements(catalog);
    const auto request = take(make_si_request(catalog, rows, admitted, plan, nullptr));
    const auto serial = take(convert_case_to_si(schema, catalog, request, admitted, plan, nullptr));
    const auto expected = serialize(serial);
    const auto worker = [&] {
        core::ResourcePlan worker_plan{budget};
        SiDiagnostic diagnostic;
        return convert_case_to_si(schema, catalog, request, admitted, worker_plan, &diagnostic);
    };
    auto first = std::async(std::launch::async, worker);
    auto second = std::async(std::launch::async, worker);
    const auto a = take(first.get());
    const auto b = take(second.get());
    CHECK(serialize(a) == expected);
    CHECK(serialize(b) == expected);
    CHECK(serialize(serial) == expected);
}

TEST_CASE("SI exact dyadic fallback propagates its arithmetic admission and allocation failures",
          "[negative][verification][IO4-01][N4-02]") {
    // 1000 * 2^-200 mm. The finite decimal coefficient has many more than
    // 256 bits before its power of ten is applied. The SI value is dyadic.
    constexpr std::string_view long_length =
        "622301527786114170714406405378012424059025216872116713310111661478969883403538344118394482"
        "31257136169569665895551224821247160434722900390625000e-200";
    const auto schema = source(long_length);
    core::ResourcePlan plan{budget};
    const auto admitted = limits(generous());
    const auto catalog = take(make_si_catalog(schema, admitted, plan, nullptr));
    const auto ordinary_rows = requirements(catalog);
    auto exact_rows = ordinary_rows;
    std::optional<SiQuantityKey> exact_key;
    for (std::uint64_t i = 0; i < take(catalog.size()); ++i) {
        const auto entry = take(catalog.entry(i));
        if (entry.field == "reference_length") {
            auto& row = exact_rows[static_cast<std::size_t>(i)];
            row.absolute = SiDecimalInput{"0"};
            row.relative = SiDecimalInput{"0"};
            REQUIRE_FALSE(exact_key);
            exact_key = entry.key;
        }
    }
    REQUIRE(exact_key);
    const auto ordinary = take(make_si_request(catalog, ordinary_rows, admitted, plan, nullptr));
    const auto exact_request = take(make_si_request(catalog, exact_rows, admitted, plan, nullptr));
    const auto allocation_count = [&](const SiRequest& request) {
        const auto result = [&] {
            const FaultScope count{-2};
            return convert_case_to_si(schema, catalog, request, admitted, plan, nullptr);
        }();
        REQUIRE(result);
        return attempts.load(std::memory_order_relaxed);
    };
    const auto ordinary_sites = allocation_count(ordinary);
    const auto exact_sites = allocation_count(exact_request);
    CHECK(exact_sites > ordinary_sites);
    const auto old =
        take(convert_case_to_si(schema, catalog, exact_request, admitted, plan, nullptr));
    const auto unchanged = serialize(old);
    const auto has_token = [](const SiExactInput& value, std::string_view expected) {
        const auto* decimal = std::get_if<SiDecimalInput>(&value);
        return decimal && decimal->token == expected;
    };
    std::size_t diagnosed_exact_allocations = 0;
    fail_allocations(
        [&](const SiLimits& resource, core::ResourcePlan& usage, SiDiagnostic* diagnostic) {
            return convert_case_to_si(schema, catalog, exact_request, resource, usage, diagnostic);
        },
        [&](const SiDiagnostic& diagnostic) {
            CHECK(serialize(old) == unchanged);
            const auto detail = diagnostic.detail();
            if (detail && detail->key == exact_key) {
                ++diagnosed_exact_allocations;
                const bool complete = detail->request_state == SiRequestState::Complete;
                CHECK(complete);
                if (const auto wanted = diagnostic.requirement()) {
                    const bool scale_matches = has_token(wanted->scale.value, "1");
                    const bool absolute_matches = has_token(wanted->absolute, "0");
                    const bool relative_matches = has_token(wanted->relative, "0");
                    CHECK(scale_matches);
                    CHECK(absolute_matches);
                    CHECK(relative_matches);
                } else
                    FAIL_CHECK("The arithmetic failure lost the original precision requirement");
            }
        });
    CHECK(diagnosed_exact_allocations > 0);

    auto small_exact = admitted.exact();
    small_exact.maximum_value_bits = 1;
    const auto short_bits = convert_case_to_si(schema, catalog, exact_request,
                                               SiLimits{generous(), small_exact}, plan, nullptr);
    REQUIRE_FALSE(short_bits);
    CHECK(short_bits.error().code == core::FailureCode::ResourceExceeded);

    auto no_external = admitted.exact();
    no_external.external_memory = structure::ExactEstimateAdmission{0};
    const auto short_external = convert_case_to_si(
        schema, catalog, exact_request, SiLimits{generous(), no_external}, plan, nullptr);
    REQUIRE_FALSE(short_external);
    CHECK(short_external.error().code == core::FailureCode::ResourceExceeded);
    CHECK(serialize(old) == unchanged);
}
