#include "../structure/exact_polynomial_detail.hpp"
#include "si_case_detail.hpp"
#include "si_case_plan.hpp"

#include <algorithm>

namespace uvlas::io::detail {
namespace {
struct SourceForecast {
    std::uint64_t quantities = 0, path = 0, depth = 0, role_work = 0;
};
SourceForecast source_forecast(SiContext& c, const InputDocumentData& source,
                               std::span<const SchemaNode> nodes) {
    c.require(SiQuota::SourceBytes, source.source.size());
    c.require(SiQuota::InputNodes, source.node_count);
    c.require(SiQuota::InputMembers, source.member_count);
    SourceForecast result;
    std::uint64_t digits = 0;
    for (std::uint64_t n = 0; n < nodes.size(); ++n) {
        c.work.step();
        c.current = n;
        const auto& node = nodes[static_cast<std::size_t>(n)];
        if (!node.assigned)
            continue;
        const auto kind = node_kind(source, n);
        std::uint64_t level = kind == InputNodeKind::Array || kind == InputNodeKind::Object ? 1 : 0,
                      path = 0, ancestor = n;
        while (nodes[static_cast<std::size_t>(ancestor)].parent) {
            c.work.step();
            const auto& info = nodes[static_cast<std::size_t>(ancestor)];
            level = c.add(level, 1);
            path = c.add(path, 1);
            if (info.parent_member) {
                const auto& range =
                    source.members.values()[static_cast<std::size_t>(*info.parent_member)].key;
                for (const char ch :
                     std::string_view{source.text.values().data() + range.offset, range.size}) {
                    c.work.step();
                    path = c.add(path, (ch == '~' || ch == '/') ? 2 : 1);
                }
            } else if (info.parent_array_index) {
                auto index = *info.parent_array_index;
                do {
                    c.work.step();
                    path = c.add(path, 1);
                    index /= 10;
                } while (index);
            }
            ancestor = *info.parent;
        }
        result.depth = std::max(result.depth, level);
        result.path = std::max(result.path, path);
        c.require(SiQuota::InputDepth, result.depth);
        if (!si_numeric_node(node))
            continue;
        result.quantities = c.add(result.quantities, 1);
        c.require(SiQuota::Quantities, result.quantities);
        const auto check = [&](const DecimalRange& number) {
            digits = c.add(digits, number.count);
            c.require(SiQuota::SourceDigits, digits);
            c.require(SiQuota::ExponentMagnitude, si::magnitude(number.exponent));
        };
        if (node_kind(source, n) == InputNodeKind::Object) {
            check(node.rational_numerator);
            check(node.rational_denominator);
        } else
            check(node.number);
        for (const auto exponent : node.dimensions.exponents)
            c.require(SiQuota::UnitPower, si::magnitude(exponent));
        const auto before = c.work.used;
        (void)si_entry(c, nodes, n);
        result.role_work = c.add(result.role_work, c.work.used - before);
    }
    c.current.reset();
    return result;
}
template <class T>
void copy(SiContext& c, const core::PlannedBuffer<T>& input, core::PlannedBuffer<T>& target) {
    target = c.allocate<T>(input.size());
    c.work.step(c.multiply(input.size(), sizeof(T)));
    std::ranges::copy(input.values(), target.values().begin());
}
template <class T, class Action>
core::Result<T> entry(SiContext& c, Action&& action) noexcept {
    try {
        return action();
    } catch (const core::Failure& f) {
        if (c.slot && !c.slot->detail()) {
            try {
                c.reject(f.code == core::FailureCode::Cancelled ? SiReason::Cancelled
                                                                : SiReason::QuotaExceeded,
                         f.code);
            } catch (...) {
                // Diagnostic failure must not replace the original value failure.
                return std::unexpected(f);
            }
        }
        return std::unexpected(f);
    } catch (const si::ArithmeticFailure& f) {
        try {
            c.reject(f.code == core::FailureCode::Cancelled ? SiReason::Cancelled
                     : f.code == core::FailureCode::ResourceExceeded
                         ? SiReason::QuotaExceeded
                         : SiReason::UnresolvedPrecision,
                     f.code);
        } catch (const core::Failure& failure) {
            return std::unexpected(failure);
        } catch (...) {
            return std::unexpected(core::Failure{f.code, core::Operation::PrepareModel});
        }
    } catch (const std::bad_alloc&) {
        try {
            c.reject(SiReason::Allocation, core::FailureCode::ResourceExceeded);
        } catch (...) {
            return std::unexpected(
                core::Failure{core::FailureCode::ResourceExceeded, core::Operation::PrepareModel});
        }
    } catch (...) {
        try {
            c.reject(SiReason::UnexpectedException);
        } catch (...) {
            return std::unexpected(
                core::Failure{core::FailureCode::InvalidInput, core::Operation::PrepareModel});
        }
    }
    return std::unexpected(
        core::Failure{core::FailureCode::InvalidInput, core::Operation::PrepareModel});
}
std::uint64_t text_bytes(SiContext& c, const SiExactInput& input) {
    if (const auto* d = std::get_if<SiDecimalInput>(&input))
        return d->token.size();
    const auto& r = std::get<SiRationalInput>(input);
    return c.add(r.numerator.size(), r.denominator.size());
}
SiExactStorage copy_exact(SiContext& c, const SiExactInput& input, std::span<char> text,
                          std::uint64_t& offset) {
    const auto save = [&](std::string_view token) {
        if (offset > text.size() || token.size() > text.size() - offset)
            c.reject(SiReason::InvalidRequirement);
        const SiTextRange range{offset, token.size()};
        c.work.step(token.size());
        std::ranges::copy(token, text.begin() + static_cast<std::ptrdiff_t>(offset));
        offset += token.size();
        return range;
    };
    if (const auto* d = std::get_if<SiDecimalInput>(&input))
        return {SiExactStorage::Decimal{save(d->token)}};
    const auto& r = std::get<SiRationalInput>(input);
    return {SiExactStorage::Rational{save(r.numerator), save(r.denominator)}};
}
} // namespace
} // namespace uvlas::io::detail

namespace uvlas::io {
core::Result<SiCatalog> make_si_catalog(const SchemaCase& input, const SiLimits& limits,
                                        core::ResourcePlan& plan, SiDiagnostic* slot) noexcept {
    detail::SiContext c{limits, plan, slot};
    return detail::entry<SiCatalog>(c, [&] {
        if (std::holds_alternative<structure::ExactRequireProvenBound>(
                limits.exact().external_memory))
            c.reject(SiReason::UnsupportedMemory, core::FailureCode::Unsupported);
        const auto* source = detail::SiFactory::storage(input);
        if (!source || !source->snapshot.size())
            c.reject(SiReason::InvalidOwner);
        const auto identity = input.identity();
        if (!identity)
            throw identity.error();
        c.source = identity->source;
        c.snapshot = &source->snapshot.values()[0];
        c.nodes = source->nodes.values();
        const auto counts = detail::source_forecast(c, *c.snapshot, c.nodes);
        auto p = c.add(c.add(sizeof(detail::SiCatalogData), sizeof(detail::InputDocumentData)),
                       sizeof(detail::SiKeyOwner));
        p = c.add(p, detail::si_snapshot_bytes(c, *c.snapshot));
        p = c.add(p, c.multiply(source->nodes.size(), sizeof(detail::SchemaNode)));
        p = c.add(p, source->digits.size());
        p = c.add(p, c.multiply(counts.quantities, sizeof(detail::SiEntryStorage)));
        const auto work = c.add(c.multiply(p, 2), c.add(source->nodes.size(), counts.role_work));
        c.admit(p, 0, work, counts.path);
        auto output = c.allocate<detail::SiCatalogData>(1);
        auto& result = output.values()[0];
        auto token = c.allocate<detail::SiKeyOwner>(1);
        const auto* token_address = token.values().data();
        result.key_owner = std::shared_ptr<const detail::SiKeyOwner>{
            token_address, [ownership = std::move(token)](
                               const detail::SiKeyOwner*) mutable noexcept { ownership = {}; }};
        result.identity = *identity;
        result.units = detail::si_units(*source);
        result.maximum_path_bytes = counts.path;
        result.snapshot = c.allocate<detail::InputDocumentData>(1);
        detail::si_clone_snapshot(c, *c.snapshot, result.snapshot.values()[0]);
        detail::copy(c, source->nodes, result.nodes);
        detail::copy(c, source->digits, result.digits);
        result.entries = c.allocate<detail::SiEntryStorage>(counts.quantities);
        std::uint64_t row = 0;
        for (std::uint64_t node = 0; node < source->nodes.size(); ++node) {
            c.work.step();
            if (detail::si_numeric_node(source->nodes.values()[static_cast<std::size_t>(node)]))
                result.entries.values()[static_cast<std::size_t>(row++)] =
                    detail::si_entry(c, c.nodes, node);
        }
        c.finish();
        result.report = c.report();
        return detail::SiFactory::catalog(std::move(output));
    });
}
core::Result<SiRequest> make_si_request(const SiCatalog& catalog,
                                        std::span<const SiRequirementRow> rows,
                                        const SiLimits& limits, core::ResourcePlan& plan,
                                        SiDiagnostic* slot) noexcept {
    detail::SiContext c{limits, plan, slot};
    return detail::entry<SiRequest>(c, [&] {
        if (std::holds_alternative<structure::ExactRequireProvenBound>(
                limits.exact().external_memory))
            c.reject(SiReason::UnsupportedMemory, core::FailureCode::Unsupported);
        const auto* source = detail::SiFactory::storage(catalog);
        if (!source || !source->identity)
            c.reject(SiReason::InvalidOwner);
        c.source = source->identity->source;
        c.snapshot = &source->snapshot.values()[0];
        c.nodes = source->nodes.values();
        const auto counts = detail::source_forecast(c, *c.snapshot, c.nodes);
        if (rows.size() != source->entries.size())
            c.reject(SiReason::MissingRequirement);
        c.require(SiQuota::Quantities, rows.size());
        std::uint64_t text = 0, tokens = 0, digits = 0, diagnostic_request = 0;
        for (const auto& row : rows) {
            c.work.step();
            c.key = row.key;
            c.requirement = SiRequirementView{row.scale, row.absolute, row.relative};
            const auto previous_text = text;
            if (!detail::same_identity(row.key.source(), *source->identity))
                c.reject(SiReason::IdentityMismatch);
            if (!detail::SiFactory::owns_key(*source, row.key))
                c.reject(SiReason::IdentityMismatch);
            const auto key_node = detail::SiFactory::node(row.key);
            if (key_node >= source->nodes.size() ||
                !detail::si_numeric_node(
                    source->nodes.values()[static_cast<std::size_t>(key_node)]))
                c.reject(SiReason::IdentityMismatch);
            c.current = key_node;
            c.dimensions = source->nodes.values()[static_cast<std::size_t>(key_node)].dimensions;
            if (!(row.scale.dimensions == *c.dimensions))
                c.reject(SiReason::DimensionMismatch);
            for (const auto* value : {&row.scale.value, &row.absolute, &row.relative}) {
                text = c.add(text, detail::text_bytes(c, *value));
                c.require(SiQuota::RequestTokenBytes, text);
                tokens = c.add(tokens, std::holds_alternative<SiDecimalInput>(*value) ? 1 : 2);
                const auto count_digits = [&](std::string_view raw) {
                    const bool negative = !raw.empty() && raw.front() == '-';
                    bool exponent = false, nonzero = false;
                    for (const char ch : raw) {
                        c.work.step();
                        if (ch == 'e' || ch == 'E')
                            exponent = true;
                        if (!exponent && ch >= '1' && ch <= '9')
                            nonzero = true;
                        if (ch >= '0' && ch <= '9') {
                            digits = c.add(digits, 1);
                            c.require(SiQuota::RequestDigits, digits);
                        }
                    }
                    return nonzero ? (negative ? -1 : 1) : 0;
                };
                int sign = 0;
                if (const auto* d = std::get_if<SiDecimalInput>(value))
                    sign = count_digits(d->token);
                else {
                    const auto& r = std::get<SiRationalInput>(*value);
                    sign = count_digits(r.numerator);
                    if (count_digits(r.denominator) <= 0)
                        c.reject(SiReason::InvalidExactInput);
                }
                // Exact signs are known from lexical digits even when the magnitude
                // is outside the interval exponent range. Never misclassify a negative budget.
                if (sign < 0 || (value == &row.scale.value && sign == 0))
                    c.reject(SiReason::InvalidRequirement);
            }
            diagnostic_request = std::max(diagnostic_request, text - previous_text);
        }
        auto p = c.add(sizeof(detail::SiRequestData), text);
        p = c.add(p, c.multiply(rows.size(), sizeof(detail::SiRequirementStorage)));
        const auto allocated_scratch =
            c.add(c.multiply(source->nodes.size(), sizeof(std::uint64_t)),
                  sizeof(detail::si::ProgramScratch));
        const auto scratch = c.add(allocated_scratch, detail::si_fixed_stack_bytes);
        // One raw token visit plus decimal multiply/add, at most 64 power bits per token.
        constexpr std::uint64_t token_power_work =
            detail::si::power_work_bound(std::numeric_limits<std::uint64_t>::max()) + 2;
        const auto work =
            c.add(c.add(c.multiply(p, 2), allocated_scratch),
                  c.add(c.multiply(text, detail::si::digit_program.work + 1),
                        c.add(c.multiply(tokens, token_power_work),
                              c.add(c.multiply(rows.size(), 9), source->nodes.size()))));
        c.admit(p, scratch, work, counts.path, allocated_scratch, 0, diagnostic_request);
        c.automatic_scratch = detail::si_fixed_stack_bytes;
        c.program = c.allocate<detail::si::ProgramScratch>(1, true);
        auto index = c.allocate<std::uint64_t>(source->nodes.size(), true);
        c.work.step(index.size());
        std::ranges::fill(index.values(), std::numeric_limits<std::uint64_t>::max());
        auto output = c.allocate<detail::SiRequestData>(1);
        auto& result = output.values()[0];
        result.identity = source->identity;
        result.key_owner = source->key_owner;
        result.rows = c.allocate<detail::SiRequirementStorage>(rows.size());
        result.text = c.allocate<char>(text);
        for (std::uint64_t row = 0; row < source->entries.size(); ++row) {
            c.work.step();
            index.values()[static_cast<std::size_t>(
                source->entries.values()[static_cast<std::size_t>(row)].node)] = row;
        }
        std::uint64_t offset = 0;
        for (const auto& row : rows) {
            c.work.step();
            c.key = row.key;
            c.requirement = SiRequirementView{row.scale, row.absolute, row.relative};
            const auto node = detail::SiFactory::node(row.key);
            c.current = node;
            if (node >= index.size())
                c.reject(SiReason::IdentityMismatch);
            const auto target = index.values()[static_cast<std::size_t>(node)];
            if (target >= source->entries.size())
                c.reject(SiReason::DuplicateRequirement);
            const auto entry = detail::si_catalog_entry(*source, target);
            if (!(row.key == entry.key))
                c.reject(SiReason::IdentityMismatch);
            c.dimensions = entry.recipe.dimensions;
            if (!(row.scale.dimensions == entry.recipe.dimensions))
                c.reject(SiReason::DimensionMismatch);
            const auto scale = detail::si_exact_interval(c, row.scale.value);
            const auto absolute = detail::si_exact_interval(c, row.absolute);
            const auto relative = detail::si_exact_interval(c, row.relative);
            c.work.step(3);
            if (detail::si::compare(scale.lower, {}) <= 0 ||
                detail::si::compare(absolute.lower, {}) < 0 ||
                detail::si::compare(relative.lower, {}) < 0)
                c.reject(SiReason::InvalidRequirement);
            result.rows.values()[static_cast<std::size_t>(target)] = {
                node, row.scale.dimensions,
                detail::copy_exact(c, row.scale.value, result.text.values(), offset),
                detail::copy_exact(c, row.absolute, result.text.values(), offset),
                detail::copy_exact(c, row.relative, result.text.values(), offset)};
            index.values()[static_cast<std::size_t>(node)] =
                std::numeric_limits<std::uint64_t>::max();
        }
        c.finish();
        result.report = c.report();
        return detail::SiFactory::request(std::move(output));
    });
}
} // namespace uvlas::io

namespace uvlas::io {
core::Result<SiCase> convert_case_to_si(const SchemaCase& input, const SiCatalog& catalog,
                                        const SiRequest& request, const SiLimits& limits,
                                        core::ResourcePlan& plan, SiDiagnostic* slot) noexcept {
    detail::SiContext c{limits, plan, slot};
    return detail::entry<SiCase>(c, [&] {
        const auto* source = detail::SiFactory::storage(input);
        const auto* inventory = detail::SiFactory::storage(catalog);
        const auto* requirements = detail::SiFactory::storage(request);
        if (!source || !inventory || !requirements || !inventory->identity ||
            !requirements->identity)
            c.reject(SiReason::InvalidOwner);
        const auto identity = input.identity();
        if (!identity)
            throw identity.error();
        c.source = identity->source;
        c.snapshot = &source->snapshot.values()[0];
        c.nodes = source->nodes.values();
        if (inventory->key_owner != requirements->key_owner ||
            !detail::same_identity(*identity, *inventory->identity) ||
            !detail::same_identity(*identity, *requirements->identity))
            c.reject(SiReason::IdentityMismatch);
        const auto counts = detail::source_forecast(c, *c.snapshot, c.nodes);
        const auto count = inventory->entries.size();
        if (counts.quantities != count || requirements->rows.size() != count)
            c.reject(SiReason::MissingRequirement);
        c.require(SiQuota::RequestTokenBytes, requirements->text.size());
        if (std::holds_alternative<structure::ExactRequireProvenBound>(
                limits.exact().external_memory))
            c.reject(SiReason::UnsupportedMemory, core::FailureCode::Unsupported);
        auto p = c.add(sizeof(detail::SiData), detail::si_schema_clone_bytes(c, *source));
        p = c.add(p, c.multiply(count, sizeof(detail::SiEntryStorage) +
                                           sizeof(detail::SiRequirementStorage) +
                                           sizeof(SiEnclosure)));
        p = c.add(p, requirements->text.size());
        p = c.add(p, c.multiply(source->nodes.size(), sizeof(std::uint64_t)));
        const auto allocated_scratch = sizeof(detail::si::ProgramScratch);
        auto work =
            c.add(c.add(c.multiply(p, 2), allocated_scratch), detail::si_rebind_work(c, *source));
        work =
            c.add(work, c.multiply(requirements->text.size(), detail::si::digit_program.work + 1));
        std::uint64_t exact_cumulative = 0, exact_peak = 0, request_tokens = 0, request_digits = 0,
                      diagnostic_request = 0;
        bool planned_pi = false;
        for (std::uint64_t row = 0; row < count; ++row) {
            c.work.step();
            const auto& entry = inventory->entries.values()[static_cast<std::size_t>(row)];
            const auto& typed = source->nodes.values()[static_cast<std::size_t>(entry.node)];
            const auto recipe = detail::si_recipe(*c.snapshot, c.nodes, source->digits.values(),
                                                  entry.node, inventory->units);
            work =
                c.add(work, detail::si_quantity_work(c, typed,
                                                     detail::node_kind(*c.snapshot, entry.node) ==
                                                         InputNodeKind::Object,
                                                     inventory->units));
            const auto exact = detail::si_exact_leaf_plan(c, recipe);
            exact_cumulative = c.add(exact_cumulative, exact.cumulative);
            exact_peak = std::max(exact_peak, exact.simultaneous);
            work = c.add(work, c.add(exact.work, 24));
            planned_pi = planned_pi || recipe.pi_power != 0;
            const auto& required = requirements->rows.values()[static_cast<std::size_t>(row)];
            const auto required_view =
                detail::si_requirement_view(required, requirements->text.values());
            diagnostic_request = std::max(
                diagnostic_request, c.add(c.add(detail::text_bytes(c, required_view.scale.value),
                                                detail::text_bytes(c, required_view.absolute)),
                                          detail::text_bytes(c, required_view.relative)));
            for (const auto* value : {&required.scale, &required.absolute, &required.relative})
                request_tokens = c.add(
                    request_tokens,
                    std::holds_alternative<detail::SiExactStorage::Decimal>(value->value) ? 1 : 2);
        }
        for (const char digit : requirements->text.values()) {
            c.work.step();
            if (digit >= '0' && digit <= '9') {
                request_digits = c.add(request_digits, 1);
                c.require(SiQuota::RequestDigits, request_digits);
            }
        }
        work = c.add(work, c.multiply(request_tokens, detail::si_token_work));
        work = c.add(work, c.multiply(count, 3)); // rational request divisions.
        work = c.add(work, source->nodes.size()); // flat index initialization.
        if (planned_pi)
            work = c.add(work, detail::si::pi_work);
        const auto scratch =
            c.add(c.add(allocated_scratch, detail::si_fixed_stack_bytes), exact_peak);
        c.admit(p, scratch, work, counts.path, allocated_scratch, exact_cumulative,
                diagnostic_request);
        c.automatic_scratch = detail::si_fixed_stack_bytes;
        auto output = c.allocate<detail::SiData>(1);
        auto& result = output.values()[0];
        result.key_owner = inventory->key_owner;
        c.program = c.allocate<detail::si::ProgramScratch>(1, true);
        result.units = inventory->units;
        result.source = detail::si_clone_schema(c, *source);
        detail::copy(c, inventory->entries, result.entries);
        detail::copy(c, requirements->rows, result.requirements);
        detail::copy(c, requirements->text, result.request_text);
        result.converted = c.allocate<SiEnclosure>(count);
        result.node_to_row = c.allocate<std::uint64_t>(source->nodes.size());
        c.work.step(result.node_to_row.size());
        std::ranges::fill(result.node_to_row.values(), std::numeric_limits<std::uint64_t>::max());
        const auto root = result.source->root();
        if (!root)
            throw root.error();
        result.defaults = {!root->provenance, !root->initial, !root->restart};
        bool needs_pi = false;
        for (const auto& record : result.entries.values()) {
            c.work.step();
            if (result.units.angle == SiAngleUnit::Degree &&
                source->nodes.values()[static_cast<std::size_t>(record.node)]
                    .dimensions.exponents[3])
                needs_pi = true;
        }
        const auto pi = needs_pi ? detail::si::programmed_pi(c.program.values()[0], c.work)
                                 : detail::si::integer(1);
        for (std::uint64_t row = 0; row < count; ++row) {
            c.work.step();
            const auto& entry = result.entries.values()[static_cast<std::size_t>(row)];
            c.current = entry.node;
            c.key = detail::si_catalog_entry(*inventory, row).key;
            const auto recipe = detail::si_recipe(*c.snapshot, c.nodes, source->digits.values(),
                                                  entry.node, result.units);
            c.dimensions = recipe.dimensions;
            const auto& wanted = result.requirements.values()[static_cast<std::size_t>(row)];
            c.requirement = detail::si_requirement_view(wanted, result.request_text.values());
            if (wanted.source_node != entry.node || !(wanted.dimensions == recipe.dimensions))
                c.reject(SiReason::IdentityMismatch);
            result.converted.values()[static_cast<std::size_t>(row)] =
                detail::si_evaluate_quantity(c, recipe, *c.requirement, pi);
            result.node_to_row.values()[static_cast<std::size_t>(entry.node)] = row;
        }
        c.require(SiQuota::CumulativeBytes, c.resources.allocated_bytes() - c.start_charged);
        c.finish();
        result.report = c.report();
        return detail::SiFactory::result(std::move(output));
    });
}
} // namespace uvlas::io
