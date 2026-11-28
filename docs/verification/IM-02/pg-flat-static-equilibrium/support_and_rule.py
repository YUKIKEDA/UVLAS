from fractions import Fraction as Q
from math import comb, factorial
import json

from positive_rational_rules import rule


def require(ok, msg):
    if not ok:
        raise ValueError(msg)


def planar(n):
    return comb(n+2,2)


def envelope(base, power):
    old = base
    for _ in range(1,power):
        old = [max(old[i]+base[k-i] for i in range(len(old))
                   if 0 <= k-i < len(base)) for k in range(len(old)+len(base)-1)]
    return old


def factored_counts():
    planes = [envelope([9,8,1], power) for power in range(1,5)]
    supports = [sum(map(planar, degree)) for degree in planes]
    f,e,p,_ = supports
    ee = 6*(5*e+planar(18)*e+e)
    g = 9*(3*p+planar(27)*f)
    stages = [18*f*f, ee, 27*f*e, g, 27*f*f, 54*f*e, g]
    return dict(planar_total_degree_envelopes_by_depth_exponent=planes,
                F_E_P_W_supports=supports,
                per_group_pairs_by_stage=stages,
                per_group_pairs=sum(stages),
                four_groups_pairs=4*sum(stages),
                active_exact_rational_slots=2*9*f+4*6*e+3*p+2*f)


def check_tensor():
    def moments(n):
        nodes, weights = rule(n)
        return tuple(sum(w*x**k for x,w in zip(nodes,weights,strict=True)) for k in range(n+1))

    radial, tangent, depth = moments(37), moments(36), moments(8)
    count = 0
    for a in range(37):
        for b in range(37-a):
            planar_rule = tangent[b]*sum((-1)**k*comb(b+1,k)*radial[a+k] for k in range(b+2))
            for c,d in enumerate(depth):
                actual = planar_rule*d
                expected = Q(factorial(a)*factorial(b), factorial(a+b+2)*(c+1))
                require(actual==expected,f'normalized triangle×depth rule for {(a,b,c)}')
                count += 1
    return count


if __name__ == '__main__':
    output=factored_counts()
    output['factored_actual_positive_Duffy_rule_moment_checks']=check_tensor()
    output['unexpanded_positive_rule_points_per_triangle_layer']=38*37*9
    print(json.dumps(output,indent=2))
