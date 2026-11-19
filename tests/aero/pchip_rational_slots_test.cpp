#include "aero/pchip_rational_slots.hpp"

#include <boost/multiprecision/cpp_int.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <future>
#include <limits>
#include <optional>
#include <vector>
namespace {
namespace p = uvlas::aero::detail::pchip;
using Big = boost::multiprecision::cpp_int;
using Rat = boost::multiprecision::cpp_rational;
template <class T>
T take(uvlas::core::Result<T> result) {
    if (!result) {
        FAIL("Unexpected failure " << static_cast<int>(result.error().code));
    }
    return std::move(*result);
}
p::RationalLimits limits() {
    return {65536, 1048576, std::numeric_limits<std::uint64_t>::max(), 10000};
}
struct Input {
    std::vector<std::uint32_t> n, d;
    int sign = 0;
    explicit Input(const Rat& value) {
        auto numerator = boost::multiprecision::numerator(value);
        sign = numerator == 0 ? 0 : numerator < 0 ? -1 : 1;
        if (numerator < 0)
            numerator = -numerator;
        boost::multiprecision::export_bits(numerator, std::back_inserter(n), 32, false);
        if (numerator == 0)
            n.clear();
        const auto denominator = boost::multiprecision::denominator(value);
        boost::multiprecision::export_bits(denominator, std::back_inserter(d), 32, false);
    }
    [[nodiscard]] p::RationalView view() const {
        return {sign, n, d};
    }
    [[nodiscard]] p::BitBound shape() const {
        return {n.empty() ? 0
                          : static_cast<std::uint64_t>(
                                (n.size() - 1) * 32 + 32 -
                                static_cast<std::size_t>(std::countl_zero(n.back()))),
                static_cast<std::uint64_t>((d.size() - 1) * 32 + 32 -
                                           static_cast<std::size_t>(std::countl_zero(d.back())))};
    }
};
Rat actual(p::Session& session, const p::Slot& slot) {
    const auto view = take(session.view(slot));
    Big n = 0, d = 0;
    for (auto i = view.numerator.size(); i > 0; --i) {
        n <<= 32;
        n += view.numerator[i - 1];
    }
    for (auto i = view.denominator.size(); i > 0; --i) {
        d <<= 32;
        d += view.denominator[i - 1];
    }
    if (view.denominator.empty()) {
        FAIL_CHECK("Missing positive denominator");
        return Rat{0};
    }
    CHECK(view.denominator.back() != 0);
    if (view.numerator.empty()) {
        CHECK(view.sign == 0);
        CHECK(d == 1);
    } else {
        CHECK(view.numerator.back() != 0);
        CHECK((view.sign == 1 || view.sign == -1));
    }
    Big ga = n, gb = d;
    while (gb != 0) {
        const Big remainder = ga % gb;
        ga = gb;
        gb = remainder;
    }
    CHECK(ga == 1);
    if (view.sign < 0)
        n = -n;
    return Rat{n} / d;
}
void compare_pair(const Rat& a, const Rat& b) {
    const Input first{a}, second{b};
    p::Counter count{limits()};
    const auto x = take(count.input(first.shape())), y = take(count.input(second.shape()));
    for (const auto op : {p::RationalOp::Add, p::RationalOp::Subtract, p::RationalOp::Multiply})
        take(count.binary(op, x, y));
    if (b != 0)
        take(count.binary(p::RationalOp::Divide, x, y));
    take(count.unary(p::RationalOp::Clone, x));
    take(count.unary(p::RationalOp::Negate, x));
    take(count.floor(x));
    REQUIRE(count.compare(x, y));
    const auto plan = take(count.finish());
    uvlas::core::ResourcePlan memory{plan.payload_bytes()};
    auto workspace = take(p::Workspace::create(plan, memory));
    auto session = take(workspace.begin());
    const auto av = take(session.input(first.view())), bv = take(session.input(second.view()));
    CHECK(actual(session, take(session.binary(p::RationalOp::Add, av, bv))) == a + b);
    CHECK(actual(session, take(session.binary(p::RationalOp::Subtract, av, bv))) == a - b);
    CHECK(actual(session, take(session.binary(p::RationalOp::Multiply, av, bv))) == a * b);
    if (b != 0)
        CHECK(actual(session, take(session.binary(p::RationalOp::Divide, av, bv))) == a / b);
    CHECK(actual(session, take(session.unary(p::RationalOp::Clone, av))) == a);
    CHECK(actual(session, take(session.unary(p::RationalOp::Negate, av))) == -a);
    Big floor = boost::multiprecision::numerator(a) / boost::multiprecision::denominator(a);
    if (a < 0 && Rat{floor} != a)
        --floor;
    CHECK(actual(session, take(session.floor(av))) == Rat{floor});
    CHECK(take(session.compare(av, bv)) == (a < b   ? std::strong_ordering::less
                                            : a > b ? std::strong_ordering::greater
                                                    : std::strong_ordering::equal));
    const auto used = take(session.usage());
    CHECK(used.steps <= plan.reserved_steps());
    CHECK(used.slots == plan.slots());
    CHECK(memory.allocated_bytes() == plan.payload_bytes());
}
} // namespace
TEST_CASE("Private PCHIP slots preserve exact signed rational arithmetic",
          "[N4-03][N4-06][verification]") {
    for (const auto a : {-7, -1, 0, 1, 7})
        for (const auto b : {-5, -1, 0, 1, 5})
            compare_pair(Rat{a} / 3, Rat{b} / 11);
    compare_pair(Rat{1} + Rat{1} / (Big{1} << 200) + Rat{1} / (Big{1} << 400), Rat{-5} / 7);
    compare_pair(Rat{std::numeric_limits<std::int64_t>::min()},
                 Rat{std::numeric_limits<std::int64_t>::max()});
}
TEST_CASE("Private PCHIP dense products and cross products match independent big integers",
          "[N4-02][N4-06][verification]") {
    for (const auto& widths : std::array<std::pair<unsigned, unsigned>, 4>{
             {{3200, 6400}, {6400, 3200}, {19000, 9000}, {16000, 16000}}})
        compare_pair(Rat{(Big{1} << widths.first) - 1}, Rat{(Big{1} << widths.second) - 3});
}
TEST_CASE("Private PCHIP preflight includes dispatch word rounding and checked overflow",
          "[N4-02][negative]") {
    p::Counter one{{1, 1, 59, 1}};
    take(one.integer(1));
    const auto plan = take(one.finish());
    CHECK(plan.integer_bits() == 1);
    CHECK(plan.words() == 1);
    CHECK(plan.slots() == 1);
    CHECK(plan.reserved_steps() == 59);
    p::Counter short_work{{1, 1, 58, 1}};
    take(short_work.integer(1));
    auto short_plan = short_work.finish();
    REQUIRE_FALSE(short_plan);
    CHECK(short_plan.error().code == uvlas::core::FailureCode::ResourceExceeded);
    p::Counter words{limits()};
    take(words.input({32, 1}));
    CHECK(take(words.finish()).words() == 2);
    p::Counter carry{limits()};
    take(carry.input({31, 1}));
    CHECK(take(carry.finish()).words() == 1);
    p::Counter no_slots{{8, 8, 10000, 0}};
    CHECK_FALSE(no_slots.integer(0));
    p::Counter no_bits{{0, 8, 10000, 1}};
    CHECK_FALSE(no_bits.integer(0));
    p::Counter narrow{{1, 1, 10000, 8}};
    const auto value = take(narrow.integer(1));
    CHECK_FALSE(narrow.binary(p::RationalOp::Multiply, value, value));
    CHECK(take(narrow.finish()).slots() ==
          1); // Failed counting never publishes a half-updated forecast.
    p::Counter huge{{UINT64_MAX, UINT64_MAX, UINT64_MAX, UINT64_MAX}};
    take(huge.input({UINT64_MAX, 1}));
    CHECK_FALSE(huge.finish());
    p::Counter invalid{limits()};
    CHECK_FALSE(invalid.input({1, 0}));
    CHECK_FALSE(invalid.unary(p::RationalOp::Count, {1, 1}));
    const p::Counter empty{{0, 0, 0, 0}};
    const auto zero = take(empty.finish());
    CHECK(zero.reserved_steps() == 0);
    uvlas::core::ResourcePlan memory{zero.payload_bytes()};
    auto w = take(p::Workspace::create(zero, memory));
    auto session = take(w.begin());
    CHECK(take(session.usage()).steps == 0);
    CHECK_FALSE(session.integer(0));
}
TEST_CASE("Private PCHIP canonical imports and failure latch preserve older slots",
          "[N4-03][N4-08][negative]") {
    p::Counter count{limits()};
    take(count.input({64, 64}));
    take(count.input({64, 64}));
    const auto plan = take(count.finish());
    const std::array<std::uint32_t, 1> one{1}, two{2}, four{4}, zero{0};
    const std::array<std::uint32_t, 2> leading{1, 0};
    const std::array<p::RationalView, 7> bad{{{1, two, four},
                                              {1, {}, one},
                                              {0, one, one},
                                              {1, zero, one},
                                              {1, one, {}},
                                              {1, one, zero},
                                              {1, leading, one}}};
    uvlas::core::ResourcePlan memory{plan.payload_bytes()};
    auto w = take(p::Workspace::create(plan, memory));
    for (const auto& input : bad) {
        auto session = take(w.begin());
        const auto baseline = take(session.input({1, one, one}));
        const auto borrowed = take(session.view(baseline));
        const auto failure = session.input(input);
        REQUIRE_FALSE(failure);
        CHECK(failure.error().code == uvlas::core::FailureCode::InvalidInput);
        CHECK(borrowed.numerator[0] == 1);
        CHECK(borrowed.denominator[0] == 1);
        const auto repeated = session.integer(1);
        REQUIRE_FALSE(repeated);
        CHECK(repeated.error() == failure.error());
    }
}
TEST_CASE("Private PCHIP leases reject stale slots and permit independent owner moves",
          "[N4-03][N4-08][negative]") {
    p::Counter count{limits()};
    take(count.integer(1));
    const auto plan = take(count.finish());
    uvlas::core::ResourcePlan memory{plan.payload_bytes()};
    auto workspace = take(p::Workspace::create(plan, memory));
    std::optional<p::Workspace> relocated;
    auto previous = [&] {
        auto session = take(workspace.begin());
        CHECK_FALSE(workspace.begin());
        auto slot = take(session.integer(1));
        relocated.emplace(std::move(workspace));
        // These queries are explicitly supported for inactive moved-from owners.
        // NOLINTBEGIN(bugprone-use-after-move)
        CHECK_FALSE(workspace.active());
        CHECK_FALSE(workspace.begin());
        // NOLINTEND(bugprone-use-after-move)
        CHECK(actual(session, slot) == 1);
        return slot;
    }();
    if (!relocated) {
        FAIL_CHECK("Missing relocated workspace");
        return;
    }
    auto& relocated_owner = relocated.value();
    auto session = take(relocated_owner.begin());
    const auto fresh = take(session.integer(1));
    const auto retained = take(session.view(fresh));
    auto stale = session.view(previous);
    REQUIRE_FALSE(stale);
    CHECK(stale.error().code == uvlas::core::FailureCode::StaleGeneration);
    CHECK(retained.numerator[0] == 1);
    CHECK(retained.denominator[0] == 1);
}
TEST_CASE("Private PCHIP cancellation stops the operation without changing previous words",
          "[N4-08][negative]") {
    p::Counter count{limits()};
    const auto value = take(count.integer(1));
    take(count.unary(p::RationalOp::Negate, value));
    const auto plan = take(count.finish());
    uvlas::core::ResourcePlan memory{plan.payload_bytes()};
    auto workspace = take(p::Workspace::create(plan, memory));
    // MSVC request_stop() is non-const; libstdc++ declares it const.
    // NOLINTNEXTLINE(misc-const-correctness)
    std::stop_source stop;
    {
        auto session = take(workspace.begin(stop.get_token()));
        const auto a = take(session.integer(1));
        const auto borrowed = take(session.view(a));
        const auto before = take(session.usage());
        stop.request_stop();
        auto failure = session.unary(p::RationalOp::Negate, a);
        REQUIRE_FALSE(failure);
        CHECK(failure.error().code == uvlas::core::FailureCode::Cancelled);
        CHECK(take(session.usage()).steps == before.steps);
        CHECK(borrowed.numerator[0] == 1);
        CHECK(borrowed.denominator[0] == 1);
    }
    CHECK_FALSE(workspace.begin(stop.get_token()));
    auto recovered = take(workspace.begin());
    CHECK(actual(recovered, take(recovered.integer(1))) == 1);
}
TEST_CASE("Private PCHIP zero division quota errors and foreign handles fail by value",
          "[N4-08][negative]") {
    p::Counter count{limits()};
    const auto a = take(count.integer(1)), b = take(count.integer(0));
    CHECK(b.numerator == 0);
    CHECK(b.denominator == 1);
    const auto division_bound = take(count.binary(p::RationalOp::Divide, a, b));
    CHECK(division_bound.denominator == 1);
    const auto plan = take(count.finish());
    uvlas::core::ResourcePlan memory{plan.payload_bytes() * 2};
    auto one = take(p::Workspace::create(plan, memory)),
         two = take(p::Workspace::create(plan, memory));
    {
        auto session = take(one.begin());
        const auto av = take(session.integer(1)), bv = take(session.integer(0));
        const auto unchanged = take(session.view(av));
        const auto unchanged_zero = take(session.view(bv));
        const auto failure = session.binary(p::RationalOp::Divide, av, bv);
        REQUIRE_FALSE(failure);
        CHECK(failure.error().code == uvlas::core::FailureCode::InvalidInput);
        CHECK(unchanged.sign == 1);
        CHECK(unchanged.numerator[0] == 1);
        CHECK(unchanged.denominator[0] == 1);
        CHECK(unchanged_zero.sign == 0);
        CHECK(unchanged_zero.numerator.empty());
        CHECK(unchanged_zero.denominator[0] == 1);
        CHECK(take(session.usage()).steps <= plan.reserved_steps());
    }
    {
        auto s = take(one.begin()), t = take(two.begin());
        const auto local = take(s.integer(1)), foreign = take(t.integer(1));
        const auto failure = s.compare(local, foreign);
        REQUIRE_FALSE(failure);
        CHECK(failure.error().code == uvlas::core::FailureCode::StaleGeneration);
    }
    {
        auto session = take(one.begin());
        take(session.integer(1));
        take(session.integer(0));
        auto failure = session.integer(1);
        REQUIRE_FALSE(failure);
        CHECK(failure.error().code == uvlas::core::FailureCode::ResourceExceeded);
    }
}
TEST_CASE("Private PCHIP independent workspace results are identical for one two and four callers",
          "[N4-03][N4-09][property]") {
    p::Counter counter{limits()};
    const auto a = take(counter.input({1, 2})), b = take(counter.integer(-7));
    const auto product = take(counter.binary(p::RationalOp::Multiply, a, b));
    take(counter.floor(product));
    const auto plan = take(counter.finish());
    struct Encoded {
        int sign = 0;
        std::vector<std::uint32_t> n, d;
        std::uint64_t steps = 0, slots = 0, payload = 0;
        bool operator==(const Encoded&) const = default;
    };
    const auto evaluate = [&]() -> std::optional<Encoded> {
        uvlas::core::ResourcePlan memory{plan.payload_bytes()};
        auto workspace = p::Workspace::create(plan, memory);
        if (!workspace)
            return std::nullopt;
        auto session = workspace->begin();
        if (!session)
            return std::nullopt;
        const std::array<std::uint32_t, 1> n{1}, d{3};
        auto first = session->input({1, n, d}), second = session->integer(-7);
        if (!first || !second)
            return std::nullopt;
        auto multiplied = session->binary(p::RationalOp::Multiply, *first, *second);
        if (!multiplied)
            return std::nullopt;
        auto floored = session->floor(*multiplied);
        if (!floored)
            return std::nullopt;
        auto view = session->view(*floored);
        auto usage = session->usage();
        if (!view || !usage)
            return std::nullopt;
        return Encoded{view->sign,
                       {view->numerator.begin(), view->numerator.end()},
                       {view->denominator.begin(), view->denominator.end()},
                       usage->steps,
                       usage->slots,
                       memory.allocated_bytes()};
    };
    const auto baseline = evaluate();
    if (!baseline) {
        FAIL_CHECK("Missing sequential output");
        return;
    }
    const auto& expected = baseline.value();
    CHECK(expected.sign == -1);
    CHECK(expected.n == std::vector<std::uint32_t>{3});
    CHECK(expected.d == std::vector<std::uint32_t>{1});
    for (const auto workers : std::array<std::size_t, 3>{1, 2, 4}) {
        std::vector<std::future<std::optional<Encoded>>> pending;
        pending.reserve(workers);
        for (std::size_t i = 0; i < workers; ++i)
            pending.push_back(std::async(std::launch::async, evaluate));
        for (auto& future : pending) {
            const auto result = future.get();
            if (!result) {
                FAIL_CHECK("Missing threaded output");
                continue;
            }
            const auto& actual_result = result.value();
            CHECK(actual_result == expected);
        }
    }
}

TEST_CASE("Private PCHIP every opcode stays within its literal limb visit forecast",
          "[N4-02][N4-06][verification]") {
    const Input source{Rat{-7} / 3};
    for (const auto op : {p::RationalOp::Input, p::RationalOp::Clone, p::RationalOp::Negate,
                          p::RationalOp::Add, p::RationalOp::Subtract, p::RationalOp::Multiply,
                          p::RationalOp::Divide, p::RationalOp::Compare, p::RationalOp::Floor}) {
        CAPTURE(static_cast<int>(op));
        p::Counter count{limits()};
        const auto shape = take(count.input(source.shape()));
        const auto second_shape = take(count.input(source.shape()));
        if (op == p::RationalOp::Input)
            take(count.input(source.shape()));
        else if (op == p::RationalOp::Compare)
            REQUIRE(count.compare(shape, second_shape));
        else if (op == p::RationalOp::Floor)
            take(count.floor(shape));
        else if (op == p::RationalOp::Clone || op == p::RationalOp::Negate)
            take(count.unary(op, shape));
        else
            take(count.binary(op, shape, second_shape));
        const auto plan = take(count.finish());
        const auto w = plan.words(), bits = plan.integer_bits();
        // Independent transcription of the adopted primitive visit bounds, plus dispatch.
        const auto multiply = w * w + 4 * w, add = 3 * w + 3;
        const auto divide = 2 * w + bits * (3 * w + 3) + 1;
        const auto normalize = (2 * bits + 3) * divide + 4 * w + 4;
        const auto budget = [&] {
            switch (op) {
            case p::RationalOp::Input:
                return normalize + 2 * w + 4;
            case p::RationalOp::Clone:
            case p::RationalOp::Negate:
                return 2 * w + 2;
            case p::RationalOp::Add:
            case p::RationalOp::Subtract:
                return 3 * multiply + add + normalize + 1;
            case p::RationalOp::Multiply:
            case p::RationalOp::Divide:
                return 2 * multiply + normalize + 1;
            case p::RationalOp::Compare:
                return 2 * multiply + add + 1;
            case p::RationalOp::Floor:
                return divide + add + 2 * w + 4;
            case p::RationalOp::Count:
                return std::uint64_t{0};
            }
            return std::uint64_t{0};
        }();
        uvlas::core::ResourcePlan memory{plan.payload_bytes()};
        auto workspace = take(p::Workspace::create(plan, memory));
        auto session = take(workspace.begin());
        const auto first = take(session.input(source.view()));
        const auto second = take(session.input(source.view()));
        const auto before = take(session.usage()).steps;
        if (op == p::RationalOp::Input)
            take(session.input(source.view()));
        else if (op == p::RationalOp::Compare)
            take(session.compare(first, second));
        else if (op == p::RationalOp::Floor)
            take(session.floor(first));
        else if (op == p::RationalOp::Clone || op == p::RationalOp::Negate)
            take(session.unary(op, first));
        else
            take(session.binary(op, first, second));
        const auto used = take(session.usage());
        CHECK(used.steps - before <= budget);
        CHECK(used.steps <= plan.reserved_steps());
        CHECK(used.slots == plan.slots());
        CHECK(memory.allocated_bytes() == plan.payload_bytes());
    }
}
