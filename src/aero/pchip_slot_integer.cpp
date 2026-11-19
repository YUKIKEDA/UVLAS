#include "pchip_slot_private.hpp"

#include <algorithm>
namespace uvlas::aero::detail::pchip::slot_detail {
// Local status propagation; neither expected refusals nor arithmetic loops throw exceptions.
bool Engine::reject(core::FailureCode code) noexcept {
    if (!data.error)
        data.error = failure(code, data.details);
    return false;
}
core::Failure Engine::error() const noexcept {
    return data.error.value_or(failure(core::FailureCode::PrecisionUnresolved, data.details));
}
bool Engine::step(std::uint64_t count) noexcept {
    if (data.stop.stop_requested())
        return reject(core::FailureCode::Cancelled);
    if (count > data.plan.reserved_steps() - data.steps)
        return reject(core::FailureCode::ResourceExceeded);
    data.steps += count;
    return true;
}
bool Engine::enter(RationalOp op) noexcept {
    if (static_cast<std::size_t>(op) >= data.remaining.size())
        return reject(core::FailureCode::InvalidInput);
    if (!(step())) {
        return false;
    }
    auto& count = data.remaining[static_cast<std::size_t>(op)];
    if (count == 0)
        return reject(core::FailureCode::ResourceExceeded);
    --count;
    return true;
}
core::Result<Slot> Engine::allocate() noexcept {
    if (data.next == data.plan.slots()) {
        reject(core::FailureCode::ResourceExceeded);
        return std::unexpected(error());
    }
    const auto index = data.next++;
    data.values.values()[static_cast<std::size_t>(index)] = {};
    const Slot result{&data, data.epoch, index};
    auto out = value(result);
    out.denominator.words[0] = 1;
    return result;
}
bool Engine::validate(const Slot& slot) noexcept {
    if (slot.owner_ != &data || slot.epoch_ != data.epoch || slot.index_ >= data.next)
        return reject(core::FailureCode::StaleGeneration);
    return true;
}
RatRef Engine::value(const Slot& slot) noexcept {
    const auto w = static_cast<std::size_t>(data.plan.words());
    auto words = data.limbs.values().subspan(static_cast<std::size_t>(slot.index_) * 2 * w, 2 * w);
    auto& m = data.values.values()[static_cast<std::size_t>(slot.index_)];
    return {{words.first(w), m.numerator}, {words.subspan(w), m.denominator}, m.sign};
}
IntRef Engine::scratch(std::size_t i) noexcept {
    const auto w = static_cast<std::size_t>(data.plan.words());
    return {
        data.limbs.values().subspan((static_cast<std::size_t>(data.plan.slots()) * 2 + i) * w, w),
        data.scratch_used[i]};
}
RationalView Engine::view(const Slot& slot) noexcept {
    auto r = value(slot);
    return {r.sign, r.numerator.read().words, r.denominator.read().words};
}
bool Engine::check_width(IntRead r) noexcept {
    return r.bits() <= data.plan.integer_bits() || reject(core::FailureCode::ResourceExceeded);
}
bool Engine::integer(IntRef out, std::uint64_t v) noexcept {
    out.used() = 0;
    while (v != 0) {
        if (!(step())) {
            return false;
        }
        if (out.used() == out.words.size())
            return reject(core::FailureCode::ResourceExceeded);
        out.words[static_cast<std::size_t>(out.used()++)] = static_cast<std::uint32_t>(v);
        v >>= 32U;
    }
    return check_width(out.read());
}
bool Engine::copy(IntRef out, IntRead in) noexcept {
    if (in.words.size() > out.words.size())
        return reject(core::FailureCode::ResourceExceeded);
    for (std::size_t i = 0; i < in.words.size(); ++i) {
        if (!(step())) {
            return false;
        }
        out.words[i] = in.words[i];
    }
    out.used() = in.words.size();
    return check_width(out.read());
}
bool Engine::compare(IntRead a, IntRead b, int& order) noexcept {
    if (!(step())) {
        return false;
    }
    order = 0;
    if (a.words.size() != b.words.size()) {
        order = a.words.size() > b.words.size() ? 1 : -1;
        return true;
    }
    for (auto i = a.words.size(); i > 0; --i) {
        if (!(step())) {
            return false;
        }
        if (a.words[i - 1] != b.words[i - 1]) {
            order = a.words[i - 1] > b.words[i - 1] ? 1 : -1;
            return true;
        }
    }
    return true;
}
bool Engine::add(IntRef out, IntRead a, IntRead b) noexcept {
    const auto size = std::max(a.words.size(), b.words.size());
    if (size > out.words.size())
        return reject(core::FailureCode::ResourceExceeded);
    std::uint64_t carry = 0, used = 0;
    for (std::size_t i = 0; i < size; ++i) {
        if (!(step())) {
            return false;
        }
        const auto v = (i < a.words.size() ? a.words[i] : 0ULL) +
                       (i < b.words.size() ? b.words[i] : 0ULL) + carry;
        out.words[i] = static_cast<std::uint32_t>(v);
        if (out.words[i] != 0)
            used = i + 1;
        carry = v >> 32U;
    }
    if (carry != 0) {
        if (!(step())) {
            return false;
        }
        if (size == out.words.size())
            return reject(core::FailureCode::ResourceExceeded);
        out.words[size] = static_cast<std::uint32_t>(carry);
        used = size + 1;
    }
    out.used() = used;
    return check_width(out.read());
}
bool Engine::subtract(IntRef out, IntRead a, IntRead b) noexcept {
    if (a.words.size() > out.words.size())
        return reject(core::FailureCode::ResourceExceeded);
    std::uint64_t borrow = 0, used = 0;
    for (std::size_t i = 0; i < a.words.size(); ++i) {
        if (!(step())) {
            return false;
        }
        const std::uint64_t first = a.words[i];
        const auto second = (i < b.words.size() ? b.words[i] : 0ULL) + borrow;
        out.words[i] = static_cast<std::uint32_t>(first - second);
        borrow = first < second ? 1 : 0;
        if (out.words[i] != 0)
            used = i + 1;
    }
    if (borrow != 0)
        return reject(core::FailureCode::PrecisionUnresolved);
    out.used() = used;
    return check_width(out.read());
}
bool Engine::multiply(IntRef out, IntRead a, IntRead b) noexcept {
    if (a.words.empty() || b.words.empty()) {
        out.used() = 0;
        return true;
    }
    const auto length = std::min(a.words.size() + b.words.size(), out.words.size());
    for (std::size_t i = 0; i < length; ++i) {
        if (!(step())) {
            return false;
        }
        out.words[i] = 0;
    }
    for (std::size_t i = 0; i < a.words.size(); ++i) {
        std::uint64_t carry = 0;
        for (std::size_t j = 0; j < b.words.size(); ++j) {
            if (!(step())) {
                return false;
            }
            if (i + j >= out.words.size())
                return reject(core::FailureCode::ResourceExceeded);
            // (2^32-1)^2 + existing limb + carry is at most 2^64-1.
            const auto v =
                static_cast<std::uint64_t>(a.words[i]) * b.words[j] + out.words[i + j] + carry;
            out.words[i + j] = static_cast<std::uint32_t>(v);
            carry = v >> 32U;
        }
        auto at = i + b.words.size();
        while (carry != 0) {
            if (!(step())) {
                return false;
            }
            if (at >= out.words.size())
                return reject(core::FailureCode::ResourceExceeded);
            const auto v = static_cast<std::uint64_t>(out.words[at]) + carry;
            out.words[at++] = static_cast<std::uint32_t>(v);
            carry = v >> 32U;
        }
    }
    auto used = length;
    while (used > 0 && out.words[used - 1] == 0) {
        if (!(step())) {
            return false;
        }
        --used;
    }
    out.used() = used;
    return check_width(out.read());
}
bool Engine::divide(IntRef quotient, IntRef remainder, IntRead a, IntRead b) noexcept {
    if (b.words.empty())
        return reject(core::FailureCode::InvalidInput);
    for (auto& word : quotient.words) {
        if (!(step())) {
            return false;
        }
        word = 0;
    }
    quotient.used() = 0;
    remainder.used() = 0;
    for (auto bit = a.bits(); bit > 0; --bit) {
        if (!(step())) {
            return false;
        }
        const auto at = bit - 1;
        std::uint64_t carry = (a.words[static_cast<std::size_t>(at / 32U)] >> (at % 32U)) & 1U;
        const auto previous = remainder.used();
        for (std::size_t i = 0; i < previous; ++i) {
            if (!(step())) {
                return false;
            }
            const auto v = (static_cast<std::uint64_t>(remainder.words[i]) << 1U) | carry;
            remainder.words[i] = static_cast<std::uint32_t>(v);
            carry = v >> 32U;
        }
        if (carry != 0) {
            if (!(step())) {
                return false;
            }
            if (previous == remainder.words.size())
                return reject(core::FailureCode::ResourceExceeded);
            remainder.words[static_cast<std::size_t>(previous)] = static_cast<std::uint32_t>(carry);
            ++remainder.used();
        }
        int order = 0;
        if (!(compare(remainder.read(), b, order))) {
            return false;
        }
        if (order >= 0) {
            if (!(subtract(remainder, remainder.read(), b))) {
                return false;
            }
            const auto word = static_cast<std::size_t>(at / 32U);
            if (word >= quotient.words.size())
                return reject(core::FailureCode::ResourceExceeded);
            if (!(step())) {
                return false;
            }
            quotient.words[word] |= std::uint32_t{1} << (at % 32U);
            quotient.used() = std::max(quotient.used(), static_cast<std::uint64_t>(word + 1));
        }
    }
    return check_width(quotient.read()) && check_width(remainder.read());
}
bool Engine::normalize(RatRef out, bool canonical) noexcept {
    if (out.denominator.used() == 0)
        return reject(core::FailureCode::InvalidInput);
    if (out.numerator.used() == 0) {
        if (canonical && !(out.denominator.used() == 1 && out.denominator.words[0] == 1))
            return reject(core::FailureCode::InvalidInput);
        out.sign = 0;
        return integer(out.denominator, 1);
    }
    auto a = scratch(0), b = scratch(1), q = scratch(2), r = scratch(3);
    if (!(copy(a, out.numerator.read()))) {
        return false;
    }
    if (!(copy(b, out.denominator.read()))) {
        return false;
    }
    while (b.used() != 0) {
        if (!(step())) {
            return false;
        }
        if (!(divide(q, r, a.read(), b.read()))) {
            return false;
        }
        std::swap(a, b);
        std::swap(b, r);
    }
    const bool one = a.used() == 1 && a.words[0] == 1;
    if (canonical)
        return one || reject(core::FailureCode::InvalidInput);
    if (!one) {
        if (!(divide(q, b, out.numerator.read(), a.read()))) {
            return false;
        }
        if (b.used() != 0)
            return reject(core::FailureCode::PrecisionUnresolved);
        if (!(copy(out.numerator, q.read()))) {
            return false;
        }
        if (!(divide(q, b, out.denominator.read(), a.read()))) {
            return false;
        }
        if (b.used() != 0)
            return reject(core::FailureCode::PrecisionUnresolved);
        if (!(copy(out.denominator, q.read()))) {
            return false;
        }
    }
    return true;
}
bool Engine::import(RatRef out, RationalView in) noexcept {
    const auto check = [&](std::span<const std::uint32_t> words, bool denominator) {
        if (words.empty())
            return !denominator || reject(core::FailureCode::InvalidInput);
        if (words.size() > data.plan.words())
            return reject(core::FailureCode::ResourceExceeded);
        if (words.back() == 0)
            return reject(core::FailureCode::InvalidInput);
        return IntRead{words}.bits() <= data.plan.limits_.input_bits ||
               reject(core::FailureCode::ResourceExceeded);
    };
    if (!(check(in.numerator, false))) {
        return false;
    }
    if (!(check(in.denominator, true))) {
        return false;
    }
    if ((in.numerator.empty() && in.sign != 0) ||
        (!in.numerator.empty() && in.sign != -1 && in.sign != 1))
        return reject(core::FailureCode::InvalidInput);
    if (!(copy(out.numerator, {in.numerator}))) {
        return false;
    }
    if (!(copy(out.denominator, {in.denominator}))) {
        return false;
    }
    out.sign = in.sign;
    return normalize(out, true);
}
bool Engine::arithmetic(RationalOp op, RatRef out, RationalView a, RationalView b) noexcept {
    if (op == RationalOp::Add || op == RationalOp::Subtract) {
        auto left = scratch(0), right = scratch(1);
        if (!(multiply(left, {a.numerator}, {b.denominator}))) {
            return false;
        }
        if (!(multiply(right, {b.numerator}, {a.denominator}))) {
            return false;
        }
        const auto second = op == RationalOp::Subtract ? -b.sign : b.sign;
        if (a.sign == second) {
            if (!(add(out.numerator, left.read(), right.read()))) {
                return false;
            }
            out.sign = a.sign;
        } else {
            int order = 0;
            if (!(compare(left.read(), right.read(), order))) {
                return false;
            }
            if (order >= 0) {
                if (!(subtract(out.numerator, left.read(), right.read()))) {
                    return false;
                }
                out.sign = a.sign;
            } else {
                if (!(subtract(out.numerator, right.read(), left.read()))) {
                    return false;
                }
                out.sign = second;
            }
        }
        if (!(multiply(out.denominator, {a.denominator}, {b.denominator}))) {
            return false;
        }
    } else if (op == RationalOp::Multiply) {
        if (!(multiply(out.numerator, {a.numerator}, {b.numerator}))) {
            return false;
        }
        if (!(multiply(out.denominator, {a.denominator}, {b.denominator}))) {
            return false;
        }
        out.sign = a.sign * b.sign;
    } else if (op == RationalOp::Divide) {
        if (b.sign == 0)
            return reject(core::FailureCode::InvalidInput);
        if (!(multiply(out.numerator, {a.numerator}, {b.denominator}))) {
            return false;
        }
        if (!(multiply(out.denominator, {a.denominator}, {b.numerator}))) {
            return false;
        }
        out.sign = a.sign * b.sign;
    } else
        return reject(core::FailureCode::InvalidInput);
    return normalize(out);
}
bool Engine::compare(RationalView a, RationalView b, std::strong_ordering& out) noexcept {
    if (a.sign != b.sign) {
        out = a.sign < b.sign ? std::strong_ordering::less : std::strong_ordering::greater;
        return true;
    }
    if (a.sign == 0) {
        out = std::strong_ordering::equal;
        return true;
    }
    auto left = scratch(0), right = scratch(1);
    if (!(multiply(left, {a.numerator}, {b.denominator}))) {
        return false;
    }
    if (!(multiply(right, {b.numerator}, {a.denominator}))) {
        return false;
    }
    int order = 0;
    if (!(compare(left.read(), right.read(), order))) {
        return false;
    }
    order *= a.sign;
    out = order < 0   ? std::strong_ordering::less
          : order > 0 ? std::strong_ordering::greater
                      : std::strong_ordering::equal;
    return true;
}
bool Engine::floor(RatRef out, RationalView in) noexcept {
    auto r = scratch(0);
    if (!(divide(out.numerator, r, {in.numerator}, {in.denominator}))) {
        return false;
    }
    if (in.sign < 0 && r.used() != 0) {
        auto one = scratch(1);
        if (!(integer(one, 1))) {
            return false;
        }
        if (!(add(out.numerator, out.numerator.read(), one.read()))) {
            return false;
        }
    }
    if (!(integer(out.denominator, 1))) {
        return false;
    }
    out.sign = out.numerator.used() == 0 ? 0 : in.sign;
    return true;
}
} // namespace uvlas::aero::detail::pchip::slot_detail
