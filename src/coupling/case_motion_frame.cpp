#include "case_motion_frame.hpp"

#include "io/si_case_detail.hpp"
#include "io/si_exact_words.hpp"
#include "structure/pg_reference_polynomials.hpp"

#include <algorithm>
#include <ranges>

namespace uvlas::coupling::detail {
namespace m = motion::m;
namespace pg = structure::detail;
namespace {
using core::FailureCode;
using pg::partition_fail;
using pg::partition_take;
using structure::ExactRational;
constexpr auto prepare_operation = core::Operation::PrepareModel;
constexpr auto check_operation = core::Operation::EvaluateStep;
constexpr std::size_t axes = 3;
constexpr unsigned interval_input_bits = 32768;
constexpr io::SchemaDimension metres{{1, 0, 0, 0}}, unitless{{0, 0, 0, 0}};

template <class F>
auto guarded(core::Operation op, core::FailureDetails details, F action) noexcept
    -> std::invoke_result_t<F> {
    const auto reframe = [&](const core::Failure& inner) {
        auto context = details;
        if (inner.details.resources)
            context.resources = inner.details.resources;
        return std::unexpected(core::Failure{inner.code, op, context});
    };
    return core::guard_entry(
        op, FailureCode::PrecisionUnresolved, details, [&]() -> std::invoke_result_t<F> {
            try {
                return action();
            } catch (const pg::PartitionAbort& failure) {
                return reframe(failure.failure);
            } catch (const m::ArithmeticFailure& failure) {
                return std::unexpected(core::Failure{failure.code, op, details});
            }
        });
}

struct FrameReader {
    FrameReader(const io::SiCase& current, const CaseMotionFrameLimits& limits,
                core::ResourcePlan& plan, const std::stop_token& stop) noexcept
        : source(current),
          arithmetic_limits{limits.exact, limits.maximum_work, 0, stop},
          arithmetic{arithmetic_limits, plan} {}

    const io::SiCase& source;
    structure::PgReferenceLimits arithmetic_limits;
    pg::ReferenceArithmetic arithmetic;

    template <class F>
    std::uint64_t scan(CaseMotionField motion, F consume) {
        const auto frame = motion.record().origin.record();
        std::uint64_t origin_node = 0;
        for (std::size_t i = 0; i < axes; ++i) {
            const auto entry =
                quantity(partition_take(frame.origin.at(static_cast<core::Index>(i))), metres);
            if (i == 0)
                origin_node = io::detail::SiFactory::node(entry.key);
            consume(true, i, entry.recipe, origin_node);
        }
        for (std::size_t i = 0; i < axes; ++i) {
            arithmetic.charge();
            const auto row = partition_take(frame.rotation.at(static_cast<core::Index>(i)));
            for (std::size_t j = 0; j < axes; ++j) {
                const auto entry =
                    quantity(partition_take(row.at(static_cast<core::Index>(j))), unitless);
                consume(false, axes * i + j, entry.recipe, origin_node);
            }
        }
        return origin_node;
    }

    [[nodiscard]] ExactRational materialize(const io::SiRecipeView& recipe) {
        const auto count = std::visit(
            [](const auto& exact) -> std::uint64_t {
                if constexpr (std::is_same_v<std::remove_cvref_t<decltype(exact)>,
                                             io::SiDecimalView>)
                    return exact.digits.size();
                else
                    return pg::partition_add(exact.numerator.digits.size(),
                                             exact.denominator.digits.size());
            },
            recipe.value);
        arithmetic.charge(pg::partition_add(1, count));
        return partition_take(
            io::materialize_rational(recipe, arithmetic.leaf_limits, arithmetic.plan));
    }

    [[nodiscard]] m::Interval enclosure(const ExactRational& value) {
        if (value.numerator_bits() > interval_input_bits ||
            value.denominator_bits() > interval_input_bits)
            partition_fail(FailureCode::ResourceExceeded);
        const auto limbs = [](std::uint64_t bits) { return (bits + 31U) / 32U; };
        arithmetic.charge(
            pg::partition_add(1, limbs(value.numerator_bits()) + limbs(value.denominator_bits())));
        if (value.sign() == 0)
            return m::integer(0);
        const auto owned = partition_take(io::detail::copy_rational_words(
            value, arithmetic.plan, arithmetic_limits.cancellation));
        const auto words = owned.view();
        const auto assemble = [&](std::span<const std::uint32_t> part) {
            aero::detail::ExactInteger integer = 0;
            for (const auto word : part | std::views::reverse) {
                arithmetic.charge();
                integer <<= 32U;
                integer += word;
            }
            return integer;
        };
        auto numerator = assemble(words.numerator);
        if (words.sign < 0)
            numerator = -numerator;
        const auto denominator = assemble(words.denominator);
        m::Work interval_work{arithmetic.scalar_operations,
                              arithmetic_limits.maximum_scalar_operations,
                              arithmetic_limits.cancellation};
        const auto result =
            m::divide(m::exact(numerator, 0), m::exact(denominator, 0), interval_work);
        arithmetic.scalar_operations = interval_work.used;
        return result;
    }

private:
    template <class Field>
    io::SiQuantityView quantity(Field field, io::SchemaDimension expected) {
        arithmetic.charge();
        auto entry = source.quantity(field);
        if (!entry && entry.error().code == FailureCode::InvalidInput)
            partition_fail(FailureCode::StaleGeneration);
        auto qualified = partition_take(std::move(entry));
        if (qualified.recipe.dimensions.exponents != expected.exponents ||
            qualified.recipe.pi_power != 0)
            partition_fail(FailureCode::InvalidInput);
        return qualified;
    }
};
} // namespace

struct CaseMotionFrame::Data {
    std::array<ExactRational, axes> origin;
    std::array<ExactRational, axes * axes> rotation;
    motion::Vector origin_interval;
    motion::Matrix rotation_interval;
};

CaseMotionFrame::CaseMotionFrame(core::PlannedBuffer<Data> data, CaseMotionFrameInfo info) noexcept
    : data_(std::move(data)), info_(info) {}

CaseMotionFrame::CaseMotionFrame(CaseMotionFrame&&) noexcept = default;
CaseMotionFrame& CaseMotionFrame::operator=(CaseMotionFrame&&) noexcept = default;
CaseMotionFrame::~CaseMotionFrame() = default;

core::Result<CaseMotionFrame>
CaseMotionFrame::prepare(const io::SiCase& source, CaseMotionField field,
                         const CaseMotionFrameLimits& limits, core::ResourcePlan& plan,
                         std::stop_token stop, core::FailureDetails details) noexcept {
    return guarded(prepare_operation, details, [&]() -> core::Result<CaseMotionFrame> {
        FrameReader reader{source, limits, plan, stop};
        reader.arithmetic.charge();
        if (limits.maximum_scratch_bytes < case_motion_frame_interval_scratch_bytes) {
            details.resources = core::ResourceUsage{core::ResourceKind::Bytes,
                                                    case_motion_frame_interval_scratch_bytes,
                                                    limits.maximum_scratch_bytes};
            return std::unexpected(
                core::Failure{FailureCode::ResourceExceeded, prepare_operation, details});
        }
        const auto identity = partition_take(partition_take(source.source_case()).get().identity());
        reader.arithmetic.admit(24);
        auto data = partition_take(plan.allocate<Data>(1, prepare_operation));
        auto& result = data.values().front();
        const auto node = reader.scan(field, [&](bool position, std::size_t i,
                                                 const io::SiRecipeView& recipe, std::uint64_t) {
            (position ? result.origin[i] : result.rotation[i]) = reader.materialize(recipe);
        });

        // The PG predicate compares all entries of Q0^T Q0 with I and checks positive
        // determinant. Exact orthogonality makes a positive determinant equal to +1.
        const auto one = reader.arithmetic.constant("1");
        pg::validate_frame(reader.arithmetic, {result.rotation, result.origin, one});
        for (std::size_t i = 0; i < result.origin.size(); ++i)
            result.origin_interval[i] = reader.enclosure(result.origin[i]);
        for (std::size_t i = 0; i < result.rotation.size(); ++i)
            result.rotation_interval[i] = reader.enclosure(result.rotation[i]);
        reader.arithmetic.cancel();
        return CaseMotionFrame{std::move(data),
                               {identity, node, reader.arithmetic.scalar_operations}};
    });
}

core::Result<void> CaseMotionFrame::validate_for(const io::SiCase& source, CaseMotionField field,
                                                 const CaseMotionFrameLimits& limits,
                                                 core::ResourcePlan& plan, std::stop_token stop,
                                                 core::FailureDetails details) const noexcept {
    return guarded(check_operation, details, [&]() -> core::Result<void> {
        FrameReader reader{source, limits, plan, stop};
        reader.arithmetic.charge();
        if (data_.size() == 0)
            partition_fail(FailureCode::StaleGeneration);
        const auto identity = partition_take(partition_take(source.source_case()).get().identity());
        if (!io::detail::same_identity(identity, info_.source))
            partition_fail(FailureCode::StaleGeneration);

        reader.arithmetic.admit(24);
        const auto& saved = data_.values().front();
        static_cast<void>(
            reader.scan(field, [&](bool position, std::size_t i, const io::SiRecipeView& recipe,
                                   std::uint64_t node) {
                if (node != info_.origin_component_node)
                    partition_fail(FailureCode::StaleGeneration);
                const auto value = reader.materialize(recipe);
                const auto& expected = position ? saved.origin[i] : saved.rotation[i];
                if (reader.arithmetic.compare(value, expected) != 0)
                    partition_fail(FailureCode::StaleGeneration);
            }));
        reader.arithmetic.cancel();
        return {};
    });
}

core::Result<void> CaseMotionFrame::match(const CaseMotionFrame& current, std::stop_token stop,
                                          core::FailureDetails details) const noexcept {
    const core::Failure cancelled{FailureCode::Cancelled, check_operation, details};
    const core::Failure stale{FailureCode::StaleGeneration, check_operation, details};
    if (stop.stop_requested())
        return std::unexpected(cancelled);
    if (data_.size() == 0 || current.data_.size() == 0 ||
        !io::detail::same_identity(info_.source, current.info_.source) ||
        info_.origin_component_node != current.info_.origin_component_node)
        return std::unexpected(stale);
    if (this == &current)
        return {};
    const auto& a = data_.values().front();
    const auto& b = current.data_.values().front();
    const auto compare = [&](const auto& left, const auto& right) -> core::Result<void> {
        for (std::size_t i = 0; i < left.size(); ++i) {
            if (stop.stop_requested())
                return std::unexpected(cancelled);
            if (left[i] != right[i])
                return std::unexpected(stale);
        }
        return {};
    };
    if (const auto origin = compare(a.origin, b.origin); !origin)
        return origin;
    if (const auto rotation = compare(a.rotation, b.rotation); !rotation)
        return rotation;
    if (stop.stop_requested())
        return std::unexpected(cancelled);
    return {};
}

core::Result<CaseMotionFrameView> CaseMotionFrame::view() const& noexcept {
    if (data_.size() == 0)
        return std::unexpected(core::Failure{FailureCode::StaleGeneration, check_operation});
    const auto& saved = data_.values().front();
    return CaseMotionFrameView{saved.origin_interval, saved.rotation_interval};
}

CaseMotionFrameSnapshot::CaseMotionFrameSnapshot(
    core::PlannedBuffer<io::detail::SiRationalWords> components, CaseMotionFrameInfo info) noexcept
    : components_(std::move(components)), info_(info) {}

core::Result<CaseMotionFrameSnapshotView> CaseMotionFrameSnapshot::view() const& noexcept {
    const auto components = components_.values();
    if (components.empty())
        return std::unexpected(core::Failure{FailureCode::InvalidInput, prepare_operation});
    CaseMotionFrameSnapshotView result{info_.source, info_.origin_component_node, {}};
    for (std::size_t i = 0; i < result.components.size(); ++i) {
        const auto words = components[i].view();
        result.components[i] = {words.sign, words.numerator, words.denominator};
    }
    return result;
}

core::Result<CaseMotionFrameSnapshot>
CaseMotionFrame::snapshot(const CaseMotionFrameLimits& limits, core::ResourcePlan& memory,
                          const std::stop_token& stop,
                          core::FailureDetails details) const noexcept {
    return guarded(prepare_operation, details, [&]() -> core::Result<CaseMotionFrameSnapshot> {
        m::Work work{0, limits.maximum_work, stop};
        work.step();
        if (data_.values().empty())
            partition_fail(FailureCode::StaleGeneration);
        const auto& current = data_.values().front();
        const auto count =
            static_cast<core::Index>(current.origin.size() + current.rotation.size());
        auto storage =
            partition_take(memory.allocate<io::detail::SiRationalWords>(count, prepare_operation));
        std::size_t i = 0;
        for (const auto parts : {std::span<const ExactRational>{current.origin},
                                 std::span<const ExactRational>{current.rotation}}) {
            for (const auto& value : parts) {
                work.step();
                const auto bits = std::max(value.numerator_bits(), value.denominator_bits());
                if (bits > limits.exact.maximum_value_bits)
                    partition_fail(FailureCode::ResourceExceeded);
                const auto n = value.numerator_bits();
                const auto d = value.denominator_bits();
                for (const auto limbs :
                     {n / 32U + (n % 32U != 0 ? 1U : 0U), d / 32U + (d % 32U != 0 ? 1U : 0U)})
                    for (std::uint64_t limb = 0; limb < limbs; ++limb)
                        work.step();
                storage.values()[i++] =
                    partition_take(io::detail::copy_rational_words(value, memory, stop));
            }
        }
        work.step();
        return CaseMotionFrameSnapshot{std::move(storage), info_};
    });
}

core::Result<CaseMotionFrame>
CaseMotionFrame::restore(const CaseMotionFrameSnapshotView& saved,
                         const CaseMotionFrameLimits& limits, core::ResourcePlan& memory,
                         const std::stop_token& stop, core::FailureDetails details) const noexcept {
    return guarded(prepare_operation, details, [&]() -> core::Result<CaseMotionFrame> {
        m::Work work{0, limits.maximum_work, stop};
        work.step();
        if (data_.values().empty() || !io::detail::same_identity(saved.source, info_.source) ||
            saved.origin_component_node != info_.origin_component_node)
            partition_fail(FailureCode::StaleGeneration);
        const auto& current = data_.values().front();
        std::size_t i = 0;
        for (const auto parts : {std::span<const ExactRational>{current.origin},
                                 std::span<const ExactRational>{current.rotation}}) {
            for (const auto& value : parts) {
                work.step();
                const auto match = [&](auto left, auto right) {
                    work.step();
                    return left == right;
                };
                const auto& from_file = saved.components[i++];
                const auto bits = std::max(value.numerator_bits(), value.denominator_bits());
                if (bits > limits.exact.maximum_value_bits)
                    partition_fail(FailureCode::ResourceExceeded);
                const auto words =
                    partition_take(io::detail::copy_rational_words(value, memory, stop));
                const auto live = words.view();
                if (from_file.sign != live.sign ||
                    !std::ranges::equal(from_file.numerator, live.numerator, match) ||
                    !std::ranges::equal(from_file.denominator, live.denominator, match))
                    partition_fail(FailureCode::StaleGeneration);
            }
        }
        auto storage = partition_take(memory.allocate<Data>(1, prepare_operation));
        auto& result = storage.values().front();
        const auto clone = [&](const auto& values, auto& output) {
            for (std::size_t n = 0; n < values.size(); ++n) {
                work.step();
                output[n] = partition_take(values[n].clone(limits.exact, memory));
            }
        };
        clone(current.origin, result.origin);
        clone(current.rotation, result.rotation);
        result.origin_interval = current.origin_interval;
        result.rotation_interval = current.rotation_interval;
        work.step();
        return CaseMotionFrame{std::move(storage), info_};
    });
}

core::Result<CaseMotionFrameInfo> CaseMotionFrame::info() const noexcept {
    if (data_.size() == 0)
        return std::unexpected(core::Failure{FailureCode::StaleGeneration, check_operation});
    return info_;
}

} // namespace uvlas::coupling::detail
