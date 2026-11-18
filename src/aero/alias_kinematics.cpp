#include "global_geometry_private.hpp"

namespace uvlas::aero::detail {
using namespace global_detail;
namespace {
std::size_t find_dof(std::span<const StoredDof> table, std::uint64_t id, Run& r) {
    std::size_t lo = 0, hi = table.size();
    while (lo < hi) {
        r.step();
        const auto mid = lo + (hi - lo) / 2;
        if (table[mid].id < id)
            lo = mid + 1;
        else
            hi = mid;
    }
    if (lo == table.size() || table[lo].id != id)
        reject(core::FailureCode::InvalidInput, core::DofId{id});
    return lo;
}
bool equal_scalar(const ParameterStorage& a, const ParameterStorage& b, Run& r) {
    return compare(scalar(restore(a), 1074, r), scalar(restore(b), 1074, r), r) == 0;
}
std::pair<std::uint16_t, std::uint16_t> pack(const ExactInteger& value, Run& r) {
    auto magnitude = value < 0 ? -value : value;
    auto& limbs = r.scratch.packed;
    limbs.fill(0);
    r.step(limbs.size());
    if (magnitude == 0)
        return {0, 0};
    const auto bits = boost::multiprecision::msb(magnitude) + 1U;
    if (bits > 8192U)
        reject(core::FailureCode::PrecisionUnsupported);
    boost::multiprecision::export_bits(magnitude, limbs.begin(), 32, false);
    const auto end = static_cast<std::size_t>((bits + 31U) / 32U);
    std::size_t first = 0;
    while (first < end && limbs[first] == 0) {
        r.step();
        ++first;
    }
    return {static_cast<std::uint16_t>(first), static_cast<std::uint16_t>(end - first)};
}
bool equal_packed(const PackedExact& p, std::span<const std::uint32_t> limbs,
                  const ExactInteger& value, Run& r) {
    const auto [first, count] = pack(value, r);
    const auto sign = value == 0 ? 0 : value < 0 ? -1 : 1;
    if (p.sign != sign || p.first != first || p.count != count)
        return false;
    for (std::size_t i = 0; i < count; ++i) {
        r.step();
        if (limbs[static_cast<std::size_t>(p.offset) + i] != r.scratch.packed[first + i])
            return false;
    }
    return true;
}
void copy_assertion(const ExactVelocityInput& input, PackedExact& output, AliasKinematicsData& data,
                    std::uint64_t& cursor, Run& r) {
    if (input.limbs.empty()) {
        if (input.sign != 0 || input.first_limb != 0)
            reject(core::FailureCode::InvalidInput);
    } else if ((input.sign != -1 && input.sign != 1) || input.first_limb >= 256 ||
               input.limbs.size() > std::size_t{256} - input.first_limb ||
               input.limbs.front() == 0 || input.limbs.back() == 0)
        reject(core::FailureCode::InvalidInput);
    output = {cursor, input.first_limb, static_cast<std::uint16_t>(input.limbs.size()), input.sign};
    for (const auto limb : input.limbs) {
        r.step();
        data.limbs.values()[static_cast<std::size_t>(cursor++)] = limb;
    }
}
ExactInteger derive(const AliasKinematicsData& data, const StoredPointKinematics& point,
                    std::size_t k, Run& r) {
    auto result = scalar(restore(point.time[k]), 3222, r); // binary64 natural grid -1074 -> -4296.
    const auto entries = data.coefficients.values().subspan(
        static_cast<std::size_t>(point.start[k]), static_cast<std::size_t>(point.length[k]));
    for (const auto& entry : entries) {
        const auto column = find_dof(data.dofs.values(), entry.column, r);
        auto coefficient = scalar(restore(entry.value), 1074, r);
        auto velocity = scalar(restore(data.dofs.values()[column].velocity), 1074, r);
        result = add(result, multiply(coefficient, velocity, r), r);
    }
    return result;
}
ExactInteger derive_input(const AliasKinematicsInput& input, const PointKinematicsInput& point,
                          std::size_t k, Run& r) {
    auto result = scalar(point.explicit_time[k], 3222, r);
    for (const auto& entry : point.jacobian[k]) {
        std::size_t lo = 0, hi = input.dofs.size();
        while (lo < hi) {
            r.step();
            const auto mid = lo + (hi - lo) / 2;
            if (input.dofs[mid].id < entry.column)
                lo = mid + 1;
            else
                hi = mid;
        }
        if (lo == input.dofs.size() || input.dofs[lo].id != entry.column)
            reject(core::FailureCode::InvalidInput, entry.column);
        result = add(
            result,
            multiply(scalar(entry.value, 1074, r), scalar(input.dofs[lo].velocity, 1074, r), r), r);
    }
    return result;
}
} // namespace
GlobalResult<AliasKinematicsState>
GlobalEngine::kinematics(const OwnedGeometry& source, const AliasKinematicsInput& input,
                         GlobalLimits limits, GlobalWorkspace& workspace, core::ResourcePlan& plan,
                         std::stop_token stop) noexcept {
    const auto context = source.context();
    return guarded(context, [&]() -> GlobalResult<AliasKinematicsState> {
        Lease const lease{workspace.busy_};
        if (!source.data_)
            reject(core::FailureCode::InvalidInput);
        const auto& g = *source.data_;
        const auto count = g.points.size();
        if (input.points.size() != count)
            reject(core::FailureCode::InvalidInput);
        Run run{limits, context, workspace.scratch_.values()[0], stop};
        // At most one point visit, three row-span visits, three assertion-span
        // visits, and one owner visit. Admit this counting phase before reading it.
        run.admit(plus(times(count, 8), 1));
        std::uint64_t coefficients = 0, assertions = 0, owners = 0;
        // Counting uses only bounded span metadata, no numerical work/allocation.
        for (const auto& point : input.points) {
            run.step();
            for (const auto row : point.jacobian) {
                run.step();
                coefficients = plus(coefficients, row.size());
            }
            if (point.velocity)
                for (const auto& v : *point.velocity) {
                    run.step();
                    assertions = plus(assertions, v.limbs.size());
                }
        }
        for (std::size_t i = 0; i < count; ++i) {
            run.step();
            if (g.points.values()[i].owner == static_cast<core::Index>(i))
                owners = plus(owners, 1);
        }
        const auto limb_ceiling = plus(times(times(owners, 3), 256), assertions);
        // Per coefficient: one exact capacity-count action, one owning output
        // action, and one optional assertion action (three products). 48 word
        // decodes also cover original validation, full alias comparison, gt,
        // basis vectors and velocities, including both words in every path.
        const auto per_term =
            plus(times(3, multiply_cost), plus(times(48, decode_cost), times(32, simple_cost)));
        const auto work = plus(
            times(plus(coefficients, plus(times(count, 3), times(input.dofs.size(), 4))), per_term),
            plus(times(limb_ceiling, 4), times(plus(count, 1), 512)));
        run.admit(work);
        metadata_admit(limits.metadata, 0);
        // Count the exact retained limb bands without allocating a DAG or data
        // array first. Successful states keep only their nonzero endpoint band,
        // never 8192 bits of zero padding per physical velocity component.
        std::uint64_t limb_capacity = assertions;
        for (std::size_t i = 1; i < input.dofs.size(); ++i) {
            run.step();
            if (input.dofs[i - 1].id >= input.dofs[i].id)
                reject(core::FailureCode::InvalidInput, input.dofs[i].id);
        }
        for (const auto& point : input.points) {
            const auto ordinal =
                lookup(g.point_lookup.values(), point.point.value(), run, point.point);
            if (g.points.values()[ordinal].owner != static_cast<core::Index>(ordinal))
                continue;
            for (std::size_t k = 0; k < 3; ++k) {
                const auto value = derive_input(input, point, k, run);
                const auto band = pack(value, run);
                limb_capacity = plus(limb_capacity, band.second);
            }
        }
        if (limb_capacity > limb_ceiling)
            reject(core::FailureCode::InvalidInput);
        auto amount = bytes<AliasKinematicsData>(1);
        amount = plus(amount, bytes<StoredDof>(input.dofs.size()));
        amount = plus(amount, bytes<StoredPointKinematics>(count));
        amount = plus(amount, bytes<StoredJacobian>(coefficients));
        amount = plus(amount, bytes<std::uint32_t>(limb_capacity));
        byte_admit(amount, plan, context);
        const auto before = plan.allocated_bytes();
        core::PlannedBuffer<AliasKinematicsData> payload;
        allocate(payload, 1, plan, context);
        auto& data = payload.values()[0];
        data.geometry.emplace(retain(source));
        data.layout_generation = input.generation.value();
        allocate(data.dofs, input.dofs.size(), plan, context);
        allocate(data.points, count, plan, context);
        allocate(data.coefficients, coefficients, plan, context);
        allocate(data.limbs, limb_capacity, plan, context);
        for (std::size_t i = 0; i < input.dofs.size(); ++i) {
            run.step();
            const auto& dof = input.dofs[i];
            if (i > 0 && input.dofs[i - 1].id >= dof.id)
                reject(core::FailureCode::InvalidInput, core::DofId{dof.id});
            auto& output = data.dofs.values()[i];
            output.id = dof.id.value();
            output.basis = dof.basis.value();
            output.units = dof.units;
            bool nonzero = false;
            for (std::size_t k = 0; k < 3; ++k) {
                output.basis_vector[k] = store(dof.basis_vector[k]);
                nonzero = scalar(dof.basis_vector[k], 1074, run) != 0 || nonzero;
            }
            if (!nonzero)
                reject(core::FailureCode::InvalidInput, core::DofId{dof.id});
            output.velocity = store(dof.velocity);
            (void)scalar(dof.velocity, 1074, run);
        }
        std::uint64_t cursor = 0, limb_cursor = 0;
        // Use row length/start only after validating point membership and uniqueness.
        // id==0 is valid; the separate seen buffer is encoded by owner=UINT64_MAX.
        for (auto& point : data.points.values())
            point.owner = std::numeric_limits<std::uint64_t>::max();
        for (const auto& point : input.points) {
            const auto index =
                lookup(g.point_lookup.values(), point.point.value(), run, point.point);
            auto& output = data.points.values()[index];
            if (output.owner != std::numeric_limits<std::uint64_t>::max())
                reject(core::FailureCode::InvalidInput, point.point);
            output.id = point.point.value();
            output.owner = static_cast<std::uint64_t>(g.points.values()[index].owner);
            output.asserted = point.velocity.has_value();
            for (std::size_t k = 0; k < 3; ++k) {
                output.start[k] = cursor;
                output.length[k] = point.jacobian[k].size();
                output.time[k] = store(point.explicit_time[k]);
                (void)scalar(point.explicit_time[k], 1074, run);
                std::optional<std::uint64_t> previous;
                for (const auto& entry : point.jacobian[k]) {
                    run.step();
                    if (previous && *previous >= entry.column.value())
                        reject(core::FailureCode::InvalidInput, point.point);
                    previous = entry.column.value();
                    (void)find_dof(data.dofs.values(), entry.column.value(), run);
                    if (scalar(entry.value, 1074, run) == 0)
                        reject(core::FailureCode::InvalidInput, point.point);
                    data.coefficients.values()[static_cast<std::size_t>(cursor++)] = {
                        entry.column.value(), store(entry.value)};
                }
                if (point.velocity)
                    copy_assertion((*point.velocity)[k], output.assertion[k], data, limb_cursor,
                                   run);
            }
        }
        // Full columns and explicit time terms precede velocity agreement; a single
        // JVP can agree for the chosen v while hiding an alias Jacobian mismatch.
        for (std::size_t i = 0; i < count; ++i) {
            const auto& p = data.points.values()[i];
            const auto& owner = data.points.values()[static_cast<std::size_t>(p.owner)];
            for (std::size_t k = 0; k < 3; ++k) {
                if (p.length[k] != owner.length[k] || !equal_scalar(p.time[k], owner.time[k], run))
                    reject(core::FailureCode::InvalidInput, PointId{p.id});
                for (std::size_t j = 0; j < p.length[k]; ++j) {
                    run.step();
                    const auto& a =
                        data.coefficients.values()[static_cast<std::size_t>(p.start[k]) + j];
                    const auto& b =
                        data.coefficients.values()[static_cast<std::size_t>(owner.start[k]) + j];
                    if (a.column != b.column || !equal_scalar(a.value, b.value, run))
                        reject(core::FailureCode::InvalidInput, PointId{p.id});
                }
            }
        }
        for (std::size_t i = 0; i < count; ++i) {
            auto& point = data.points.values()[i];
            if (point.owner != i)
                continue;
            for (std::size_t k = 0; k < 3; ++k) {
                auto value = derive(data, point, k, run);
                const auto [first, length] = pack(value, run);
                point.velocity[k] = {limb_cursor, first, length,
                                     static_cast<std::int8_t>(value == 0  ? 0
                                                              : value < 0 ? -1
                                                                          : 1)};
                for (std::size_t j = 0; j < length; ++j) {
                    run.step();
                    data.limbs.values()[static_cast<std::size_t>(limb_cursor++)] =
                        run.scratch.packed[first + j];
                }
            }
        }
        for (std::size_t i = 0; i < count; ++i) {
            auto& point = data.points.values()[i];
            const auto& owner = data.points.values()[static_cast<std::size_t>(point.owner)];
            if (point.owner != i)
                point.velocity = owner.velocity;
            if (point.asserted)
                for (std::size_t k = 0; k < 3; ++k) {
                    const auto value = derive(data, point, k, run);
                    if (!equal_packed(point.assertion[k], data.limbs.values(), value, run))
                        reject(core::FailureCode::InvalidInput, PointId{point.id});
                }
        }
        if (limb_cursor != limb_capacity)
            reject(core::FailureCode::InvalidInput);
        data.plan = {plan.allocated_bytes() - before, run.reserved, run.used, 0, 0};
        return AliasKinematicsState{std::move(payload), context};
    });
}
AliasKinematicsState::AliasKinematicsState(core::PlannedBuffer<AliasKinematicsData> p,
                                           GeometryContext c) noexcept
    : data_(std::move(p)), context_(c) {}
AliasKinematicsState::~AliasKinematicsState() = default;
AliasKinematicsState::AliasKinematicsState(AliasKinematicsState&&) noexcept = default;
AliasKinematicsState& AliasKinematicsState::operator=(AliasKinematicsState&&) noexcept = default;
bool AliasKinematicsState::active() const noexcept {
    return !data_.values().empty();
}
GlobalResult<ExactVelocityView> AliasKinematicsState::velocity(PointId id,
                                                               std::size_t k) const& noexcept {
    if (!active() || k >= 3)
        return std::unexpected(failure(core::FailureCode::InvalidInput, context_, id));
    const auto& data = data_.values()[0];
    const auto& g = *GlobalEngine::source(*data.geometry);
    const auto table = g.point_lookup.values();
    const auto found =
        std::lower_bound(table.begin(), table.end(), id.value(),
                         [](const auto& row, std::uint64_t key) { return row.id < key; });
    if (found == table.end() || found->id != id.value())
        return std::unexpected(failure(core::FailureCode::InvalidInput, context_, id));
    const auto& p = data.points.values()[static_cast<std::size_t>(found->ordinal)].velocity[k];
    return ExactVelocityView{
        p.sign, p.first, data.limbs.values().subspan(static_cast<std::size_t>(p.offset), p.count)};
}
GlobalResult<std::reference_wrapper<const OwnedGeometry>>
AliasKinematicsState::geometry() const& noexcept {
    if (!active())
        return std::unexpected(failure(core::FailureCode::InvalidInput, context_));
    return std::cref(*data_.values()[0].geometry);
}
GlobalResult<GlobalWorkPlan> AliasKinematicsState::work_plan() const noexcept {
    if (!active())
        return std::unexpected(failure(core::FailureCode::InvalidInput, context_));
    return data_.values()[0].plan;
}
} // namespace uvlas::aero::detail
