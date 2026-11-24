#pragma once

#include "nastran_support_cards.hpp"

#include <algorithm>

namespace uvlas::io::detail {

// Traverse in original term order. SPCADD takes precedence over direct terms with the same SID.
template <class Step, class Consume>
void visit_nastran_support_terms(const NastranSupportView& supports, NastranSupportSetId selected,
                                 NastranGridId requested, Step step, Consume consume) {
    const auto member_of = [&](const auto& matches) {
        return std::ranges::any_of(supports.union_members, [&](const auto& member) {
            step();
            return member.combination == selected && matches(member);
        });
    };
    const bool combination = member_of([](const auto&) { return true; });
    for (const auto& term : supports.terms) {
        step();
        const bool in_set =
            combination ? member_of([&](const auto& member) { return member.member == term.set; })
                        : term.set == selected;
        if (!in_set) {
            continue;
        }
        const bool includes = std::visit(
            [&](const auto& target) {
                if constexpr (requires { target.grid; }) {
                    return target.grid == requested;
                } else {
                    return requested >= target.first && requested <= target.last;
                }
            },
            term.target);
        if (includes) {
            consume(term);
        }
    }
}

} // namespace uvlas::io::detail
