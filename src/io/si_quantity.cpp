#include "si_case_detail.hpp"

namespace uvlas::io::detail {
// Called only after exact source/request validation and whole-call admission.
// A scalar result never creates or substitutes for public SiCase qualification.
SiEnclosure si_evaluate_quantity(SiContext& c, const SiRecipeView& recipe,
                                 const SiRequirementView& wanted, const si::Interval& pi) {
    c.dimensions = recipe.dimensions;
    c.requirement = wanted;
    si::Interval original;
    if (const auto* decimal = std::get_if<SiDecimalView>(&recipe.value))
        original = si::programmed_decimal(*decimal, c.program.values()[0], c.work);
    else {
        const auto& r = std::get<SiRationalView>(recipe.value);
        original = si::divide(si::programmed_decimal(r.numerator, c.program.values()[0], c.work),
                              si::programmed_decimal(r.denominator, c.program.values()[0], c.work),
                              c.work);
    }
    const auto truth = si::multiply(original, si_unit_interval(c, recipe, pi), c.work);
    const auto scale = si_exact_interval(c, wanted.scale.value);
    const auto absolute = si_exact_interval(c, wanted.absolute);
    const auto relative = si_exact_interval(c, wanted.relative);
    const auto budget = si::add(si::multiply(scale, absolute, c.work),
                                si::multiply(relative, si::absolute(truth), c.work), c.work);
    auto result = si::center(truth, false, c.work);
    c.achieved = result.value;
    c.work.step();
    if (si::compare(result.error, budget.lower) > 0) {
        result = si::center(truth, true, c.work);
        c.achieved = result.value;
        c.work.step();
        if (si::compare(result.error, budget.lower) > 0) {
            const auto exact = si_exact_zero_center(c, recipe, true);
            if (!exact)
                c.reject(SiReason::UnresolvedPrecision, core::FailureCode::PrecisionUnresolved);
            result.value = *exact;
            c.achieved = *exact;
        }
    }
    return result.value;
}
} // namespace uvlas::io::detail
