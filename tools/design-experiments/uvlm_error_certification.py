"""Numerical-study acceptance only; never asserts physical validation."""
import math
def estimate(values,ratio=2.,p_min=.5,p_max=2.,order_drift=.25,safety=1.5):
    if len(values)!=4 or ratio<=1 or safety<1:raise ValueError('invalid_study')
    q=[complex(x) for x in values]
    if any(not(math.isfinite(x.real) and math.isfinite(x.imag)) for x in q):return dict(resolved=False,reason='nonfinite_observation')
    d=[q[i+1]-q[i] for i in range(3)]
    if min(abs(x) for x in d)==0:return dict(resolved=False,reason='zero_difference_is_not_a_bound')
    directions=[(d[i].conjugate()*d[i+1]).real/(abs(d[i])*abs(d[i+1])) for i in range(2)]
    if min(directions)<.99:return dict(resolved=False,reason='nonmonotone_or_rotating_increment')
    orders=[math.log(abs(d[i])/abs(d[i+1]))/math.log(ratio) for i in range(2)]
    if min(orders)<p_min or max(orders)>p_max or abs(orders[0]-orders[1])>order_drift:return dict(resolved=False,reason='unqualified_asymptotic_order',orders=orders)
    p=min(orders);correction=d[-1]/(ratio**p-1);limit=q[-1]+correction;older=q[-2]+d[-2]/(ratio**orders[0]-1);roundoff=32*2.220446049250313e-16*max(abs(x) for x in q)
    return dict(resolved=True,orders=orders,finite_error_indicator=safety*abs(correction)+roundoff,limit=limit,reference_uncertainty_indicator=safety*abs(limit-older)+roundoff,rigorous_bound=False)
def certify(plan,estimates):
    if plan['scale']<=0 or plan['relative_tolerance']<0:raise ValueError('invalid_quantity_scale')
    needed=set(plan['required']);covered=set();total=0.;reasons=[]
    for e in estimates:
        groups=set(e['covers'])
        if covered&groups:reasons.append('overlapping_error_groups')
        covered|=groups
        if e.get('case_id')!=plan['case_id']:reasons.append('physical_case_mismatch')
        if not e.get('resolved',False):reasons.append('unresolved_error_group');continue
        value=e['absolute_indicator']
        if not math.isfinite(value) or value<0:raise ValueError('invalid_error_indicator')
        total+=value
        for group,cap in plan.get('absolute_group_caps',{}).items():
            if group in groups:
                if len(groups)!=1:reasons.append('individual_cap_not_separately_resolved')
                elif value>cap:reasons.append('individual_cap_exceeded')
    if covered!=needed:reasons.append('missing_or_unexpected_error_groups')
    if total>plan['scale']*plan['relative_tolerance']:reasons.append('total_budget_exceeded')
    return dict(qualified=not reasons,reasons=sorted(set(reasons)),absolute_indicator=total,relative_indicator=total/plan['scale'],physical_validation=False)
