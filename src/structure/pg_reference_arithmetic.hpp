#pragma once

#include "pg_partition_common.hpp"

#include <array>
#include <compare>
#include <string_view>
#include <utility>
#include <uvlas/structure/pg_reference_geometry.hpp>

namespace uvlas::structure::detail {

/// Fixed live slots of the largest point path (including input/output and nested
/// jet, matrix, and radical expressions). Each exact op charges its own allocation.
inline constexpr std::uint64_t reference_live_rationals = 256;

struct ReferenceArithmetic {
    const PgReferenceLimits& limits;
    core::ResourcePlan& plan;
    std::uint64_t scalar_operations = 0;
    std::uint64_t partition_visits = 0;
    ExactLimits leaf_limits = limits.algebra;
    std::uint64_t initial_bytes = plan.allocated_bytes();

    [[nodiscard]] PgReferenceCallInfo info() const noexcept {
        return {scalar_operations, partition_visits, plan.allocated_bytes() - initial_bytes};
    }

    void visit(std::uint64_t count) {
        cancel();
        if (count > limits.maximum_partition_visits - partition_visits)
            partition_fail(core::FailureCode::ResourceExceeded);
        partition_visits += count;
    }

    void cancel() const {
        if (limits.cancellation.stop_requested())
            partition_fail(core::FailureCode::Cancelled);
    }
    void charge() {
        cancel();
        if (scalar_operations == limits.maximum_scalar_operations)
            partition_fail(core::FailureCode::ResourceExceeded);
        ++scalar_operations;
    }
    void admit() {
        cancel();
        if (!std::holds_alternative<ExactEstimateAdmission>(limits.algebra.external_memory))
            partition_fail(core::FailureCode::Unsupported);
        constexpr std::uint64_t limb_bits = IntegerBackend::limb_bits;
        constexpr auto limb_bytes = sizeof(boost::multiprecision::limb_type);
        const auto limbs = limits.algebra.maximum_work_bits / limb_bits;
        if (limbs < 3)
            partition_fail(core::FailureCode::ResourceExceeded);
        const auto grown = partition_mul(4, limbs);
        // Two integers per live rational and one spare backend slot for add/sub.
        const auto extra =
            partition_mul(partition_mul(2 * reference_live_rationals + 1, grown), limb_bytes);
        partition_take(exact_rational_admission(limits.algebra, (limbs - 2) * limb_bits, extra));
        // Leave the same held-slot allowance in place when a leaf measures the
        // actual borrowed limbs of the participating operands.
        std::get<ExactEstimateAdmission>(leaf_limits.external_memory).bytes -= extra;
    }
    [[nodiscard]] ExactRational constant(std::string_view value) {
        charge();
        return partition_take(ExactRational::from_integer_pair(value, "1", leaf_limits, plan));
    }
    [[nodiscard]] ExactRational copy(const ExactRational& value) {
        charge();
        return partition_take(value.clone(leaf_limits, plan));
    }
    [[nodiscard]] ExactRational add(const ExactRational& a, const ExactRational& b) {
        charge();
        return partition_take(a.add(b, leaf_limits, plan));
    }
    [[nodiscard]] ExactRational sub(const ExactRational& a, const ExactRational& b) {
        charge();
        return partition_take(a.subtract(b, leaf_limits, plan));
    }
    [[nodiscard]] ExactRational mul(const ExactRational& a, const ExactRational& b) {
        charge();
        return partition_take(a.multiply(b, leaf_limits, plan));
    }
    [[nodiscard]] ExactRational div(const ExactRational& a, const ExactRational& b) {
        charge();
        return partition_take(a.divide(b, leaf_limits, plan));
    }
    [[nodiscard]] std::strong_ordering compare(const ExactRational& a, const ExactRational& b) {
        charge();
        return partition_take(a.compare(b, leaf_limits, plan));
    }
    [[nodiscard]] ExactRational neg(const ExactRational& a) {
        return sub({}, a);
    }
    [[nodiscard]] ExactRational abs(const ExactRational& a) {
        return a.sign() < 0 ? neg(a) : copy(a);
    }
    [[nodiscard]] ExactRational square(const ExactRational& a) {
        return mul(a, a);
    }
};

/// Arithmetic in the represented algebra a+b sqrt(g); no division in the
/// quadratic algebra is used, so a rational-square g poses no ambiguity.
struct ReferenceRadicals {
    ReferenceArithmetic& math;
    const ExactRational& g;

    [[nodiscard]] PgReferenceScalar add(const PgReferenceScalar& a, const PgReferenceScalar& b) {
        return {math.add(a.rational, b.rational), math.add(a.radical, b.radical)};
    }
    [[nodiscard]] PgReferenceScalar sub(const PgReferenceScalar& a, const PgReferenceScalar& b) {
        return {math.sub(a.rational, b.rational), math.sub(a.radical, b.radical)};
    }
    [[nodiscard]] PgReferenceScalar scale(const PgReferenceScalar& a, const ExactRational& x) {
        return {math.mul(a.rational, x), math.mul(a.radical, x)};
    }
    [[nodiscard]] PgReferenceScalar mul(const PgReferenceScalar& a, const PgReferenceScalar& b) {
        return {
            math.add(math.mul(a.rational, b.rational), math.mul(g, math.mul(a.radical, b.radical))),
            math.add(math.mul(a.rational, b.radical), math.mul(a.radical, b.rational))};
    }
    [[nodiscard]] int sign(const PgReferenceScalar& a) {
        const int p = a.rational.sign(), q = a.radical.sign();
        if (p == 0)
            return q;
        if (q == 0 || p == q)
            return p;
        const auto ratio =
            math.compare(math.square(a.rational), math.mul(g, math.square(a.radical)));
        return ratio == 0 ? 0 : (ratio > 0 ? p : q);
    }
    [[nodiscard]] PgReferenceVector cross(const PgReferenceVector& a, const PgReferenceVector& b) {
        PgReferenceVector result;
        for (std::size_t i = 0; i < result.size(); ++i) {
            const auto j = (i + 1) % result.size(), k = (i + 2) % result.size();
            result[i] = sub(mul(a[j], b[k]), mul(a[k], b[j]));
        }
        return result;
    }
    [[nodiscard]] PgReferenceScalar dot(const PgReferenceVector& a, const PgReferenceVector& b) {
        PgReferenceScalar result;
        for (std::size_t i = 0; i < a.size(); ++i)
            result = add(result, mul(a[i], b[i]));
        return result;
    }
};

} // namespace uvlas::structure::detail
