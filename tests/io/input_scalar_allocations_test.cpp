#include "../support/allocation_audit.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <new>
#include <uvlas/io/input_scalar.hpp>

namespace {
using namespace uvlas;
using namespace uvlas::io;
using namespace uvlas::test_support;
using Diagnostic = std::optional<ScalarConversionDiagnostic>;

template <class Id>
[[nodiscard]] bool exercise_id() {
    Diagnostic diagnostic;
    const auto maximum = read_physical_id<Id>("18446744073709551615", diagnostic);
    const bool maximum_ok = maximum.has_value() && maximum->value() == 18446744073709551615ULL &&
                            !diagnostic.has_value();
    const auto rounded = read_physical_id<Id>("9007199254740993", diagnostic);
    const bool rounded_ok =
        rounded.has_value() && rounded->value() == 9007199254740993ULL && !diagnostic.has_value();
    for (const auto text :
         {"", "+1", "-1", "1 ", "1.0", "1e0", "18446744073709551616", "\xef\xbc\x91"}) {
        const auto failed = read_physical_id<Id>(text, diagnostic);
        if (failed.has_value() || !diagnostic.has_value() ||
            failed.error().code != core::FailureCode::InvalidInput)
            return false;
    }
    const auto zero = read_physical_id<Id>("0000", diagnostic);
    return maximum_ok && rounded_ok && zero.has_value() && zero->value() == 0 &&
           !diagnostic.has_value();
}
} // namespace

TEST_CASE("Every strict scalar path performs no C++ allocation",
          "[IO4-01][N4-02][N4-09][property][negative]") {
    // This product path has no C heap, Eigen or external-backend operations.
    // Catch assertions and thread creation stay outside the audited interval.
    AllocationAudit audit{cpp_backend, AuditThreadDomain::AllThreads};
    const std::array domains{exercise_id<core::DofId>(),        exercise_id<core::PanelId>(),
                             exercise_id<core::WakeVertexId>(), exercise_id<core::WakeFaceId>(),
                             exercise_id<core::LayerId>(),      exercise_id<core::ObservationId>(),
                             exercise_id<core::ScopeId>(),      exercise_id<core::ComponentId>()};
    Diagnostic diagnostic;
    const auto signed_ok = read_index(NativeInteger{std::int64_t{0}}, diagnostic);
    const auto unsigned_ok =
        read_index(NativeInteger{std::uint64_t{9223372036854775807ULL}}, diagnostic);
    const auto negative = read_index(NativeInteger{std::int64_t{-1}}, diagnostic);
    const auto overflow =
        read_index(NativeInteger{std::uint64_t{18446744073709551615ULL}}, diagnostic);
    const auto allocations = audit.finish();
    for (const bool good : domains)
        CHECK(good);
    CHECK(signed_ok.has_value());
    CHECK(unsigned_ok.has_value());
    CHECK_FALSE(negative.has_value());
    CHECK_FALSE(overflow.has_value());
    CHECK(allocations.status == AllocationAuditStatus::Clean);
    CHECK(allocations.cpp_allocation_attempts == 0);
}

TEST_CASE("Scalar allocation audit observes a deliberately allocated control",
          "[N4-09][property][negative]") {
    AllocationAudit audit{cpp_backend, AuditThreadDomain::AllThreads};
    // Direct allocation calls cannot be optimized away like unused new-expressions.
    void* pointer = ::operator new(16);
    ::operator delete(pointer);
    const auto observed = audit.finish();
    CHECK(observed.status == AllocationAuditStatus::AllocationObserved);
    CHECK(observed.cpp_allocation_attempts > 0);
}
