#include "input_document_detail.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <limits>
#include <new>
#include <optional>
#include <string_view>

namespace {
// Process-local injection belongs only to this dedicated test executable.
enum class InjectedFailure : std::uint8_t { Allocation, Unexpected, Nonstandard };
struct NonstandardException {};
constinit std::atomic<std::int64_t> remaining{-1};
constinit std::atomic<std::uint64_t> attempts{0};
constinit InjectedFailure injected = InjectedFailure::Allocation;
class ForeignException final : public std::exception {
public:
    [[nodiscard]] const char* what() const noexcept override {
        return "fixture exception";
    }
};
void* allocate(std::size_t count) {
    const auto next = remaining.load(std::memory_order_relaxed);
    if (next != -1) {
        attempts.fetch_add(1, std::memory_order_relaxed);
        if (next == 0) {
            // Disable before exception construction and stack unwinding.
            remaining.store(-1, std::memory_order_relaxed);
            if (injected == InjectedFailure::Unexpected)
                throw ForeignException{};
            if (injected == InjectedFailure::Nonstandard)
                throw NonstandardException{};
            throw std::bad_alloc{};
        }
        if (next > 0)
            remaining.fetch_sub(1, std::memory_order_relaxed);
    }
    // NOLINTNEXTLINE(misc-const-correctness): operator new must return writable pointee storage.
    void* const value = std::malloc(count == 0 ? 1 : count);
    if (value == nullptr)
        throw std::bad_alloc{};
    return value;
}
struct FaultScope {
    explicit FaultScope(std::int64_t count,
                        InjectedFailure kind = InjectedFailure::Allocation) noexcept {
        attempts.store(0, std::memory_order_relaxed);
        injected = kind;
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

namespace {
using namespace uvlas;
using namespace uvlas::io;
constexpr std::uint64_t budget = 64ULL * 1024 * 1024;
constexpr std::string_view source =
    R"({"schema":"uvlas.case/1","units":null,"model":["日本語",1.0,{"x":null}],"flow":null,"analysis":null,"observations":null,"numerics":null,"resources":null,"output":null})";
InputLimits limits() {
    return {100000,
            64,
            10000,
            10000,
            10000,
            10000,
            100000,
            {ParserAllocationMode::EstimateAdmission, std::numeric_limits<std::uint64_t>::max()}};
}
core::Result<InputDocument> read(const InputLimits& request, core::ResourcePlan& resources,
                                 std::optional<InputDiagnostic>& diagnostic,
                                 core::FailureDetails details = {}) {
    return read_input_document(std::as_bytes(std::span{source.data(), source.size()}), SourceId{7},
                               request, resources, diagnostic, details);
}
} // namespace

TEST_CASE("Document allocation failures preserve input old owners and cumulative plan charges",
          "[negative][IO4-01][N4-02][N4-03][N4-08]") {
    core::ResourcePlan original_plan{budget};
    std::optional<InputDiagnostic> diagnostic;
    const auto old = read(limits(), original_plan, diagnostic);
    REQUIRE(old);
    const auto old_view = old->view();
    REQUIRE(old_view);
    const auto identity = old_view->identity;
    const auto original_bytes = std::as_bytes(std::span{source.data(), source.size()});
    core::ResourcePlan control_plan{budget};
    const auto control = [&] {
        const FaultScope count{-2};
        return read(limits(), control_plan, diagnostic);
    }();
    REQUIRE(control);
    const auto count = attempts.load(std::memory_order_relaxed);
    REQUIRE(count > 8); // Product arrays plus the actual external SAX path.
    const core::FailureDetails context{core::ModelGeneration{1},
                                       core::AcceptedStateId{2},
                                       core::TrialId{3},
                                       core::NumericalGeneration{4},
                                       core::ComponentId{5},
                                       core::ProblemId{core::PanelId{6}},
                                       std::nullopt,
                                       core::AcceptedStateId{9}};
    bool retained_successful_charge = false;
    for (std::uint64_t site = 0; site < count; ++site) {
        core::ResourcePlan plan{budget};
        const auto prior = plan.allocate<std::byte>(3, core::Operation::PrepareModel);
        REQUIRE(prior);
        const auto failed = [&] {
            const FaultScope fault{static_cast<std::int64_t>(site)};
            return read(limits(), plan, diagnostic, context);
        }();
        REQUIRE_FALSE(failed);
        CHECK(failed.error().code == core::FailureCode::ResourceExceeded);
        CHECK(failed.error().operation == core::Operation::PrepareModel);
        CHECK(failed.error().details.model == context.model);
        CHECK(failed.error().details.parent == context.parent);
        CHECK(failed.error().details.trial == context.trial);
        CHECK(failed.error().details.numerical_generation == context.numerical_generation);
        CHECK(failed.error().details.component == context.component);
        CHECK(failed.error().details.problem == context.problem);
        CHECK(failed.error().details.last_accepted == context.last_accepted);
        if (!diagnostic) {
            FAIL_CHECK("Missing allocation diagnostic");
            return;
        }
        CHECK(diagnostic->reason == InputDiagnosticReason::Allocation);
        CHECK(diagnostic->source.value() == 7);
        CHECK(plan.allocated_bytes() >= 3);
        retained_successful_charge = retained_successful_charge || plan.allocated_bytes() > 3;
        CHECK(old->view()->identity == identity);
        CHECK(std::equal(old_view->source_bytes.begin(), old_view->source_bytes.end(),
                         original_bytes.begin()));
    }
    CHECK(retained_successful_charge);
    core::ResourcePlan recovery{budget};
    CHECK(read(limits(), recovery, diagnostic));
    CHECK_FALSE(diagnostic);
    std::printf("input-document actual C++new fault sites=%llu\n",
                static_cast<unsigned long long>(count));
}

TEST_CASE(
    "Document entry contains foreign exceptions and refuses known policies without allocation",
    "[negative][IO4-01][N4-02][N4-08]") {
    core::ResourcePlan plan{budget};
    std::optional<InputDiagnostic> diagnostic;
    for (const auto kind : {InjectedFailure::Unexpected, InjectedFailure::Nonstandard}) {
        const auto foreign = [&] {
            const FaultScope fault{0, kind};
            return read(limits(), plan, diagnostic);
        }();
        REQUIRE_FALSE(foreign);
        CHECK(foreign.error().code == core::FailureCode::InvalidInput);
        if (!diagnostic) {
            FAIL_CHECK("Missing foreign-exception diagnostic");
            return;
        }
        CHECK(diagnostic->reason == InputDiagnosticReason::UnexpectedException);
        CHECK(plan.allocated_bytes() == 0);
    }
    auto strict = limits();
    strict.external_memory.mode = ParserAllocationMode::RequireProvenBound;
    const auto rejected = [&] {
        const FaultScope count{-2};
        return read(strict, plan, diagnostic);
    }();
    REQUIRE_FALSE(rejected);
    CHECK(rejected.error().code == core::FailureCode::Unsupported);
    CHECK(attempts.load(std::memory_order_relaxed) == 0);
    CHECK(plan.allocated_bytes() == 0);
    auto shortage = limits();
    shortage.external_memory.maximum_bytes = 0;
    const auto too_small = [&] {
        const FaultScope count{-2};
        return read(shortage, plan, diagnostic);
    }();
    REQUIRE_FALSE(too_small);
    CHECK(too_small.error().code == core::FailureCode::ResourceExceeded);
    CHECK(attempts.load(std::memory_order_relaxed) == 0);
}

TEST_CASE(
    "Parser message fallback classification owns no third party message and allocates nothing",
    "[unit][negative][IO4-01][N4-08]") {
    using Kind = io::detail::ParserFailureClass;
    using Code = core::FailureCode;
    // This reaches the real private classifier. It does NOT inject CRT malloc.
    for (const auto kind : {Kind::Syntax, Kind::NumberRange}) {
        const auto result = [&] {
            const FaultScope count{-2};
            return io::detail::classify_parser_failure(kind, kind == Kind::Syntax ? 101 : 406,
                                                       "Unknown exception");
        }();
        CHECK(result.code == Code::ResourceExceeded);
        CHECK(result.reason == InputDiagnosticReason::Allocation);
        CHECK(attempts.load(std::memory_order_relaxed) == 0);
    }
    const auto foreign = io::detail::classify_parser_failure(Kind::Other, 101, "Unknown exception");
    CHECK(foreign.code == Code::InvalidInput);
    CHECK(foreign.reason == InputDiagnosticReason::UnexpectedException);
    const auto range = io::detail::classify_parser_failure(
        Kind::NumberRange, 406, "[json.exception.out_of_range.406] fixed");
    CHECK(range.code == Code::Nonfinite);
    CHECK(range.reason == InputDiagnosticReason::NumberRange);
    const auto malformed =
        io::detail::classify_parser_failure(Kind::Syntax, 101, "unexpected message");
    CHECK(malformed.reason == InputDiagnosticReason::UnexpectedException);
}
