"""GC cell-break variant: Q2 midpoint nodes are not polynomial breaks.

The curvature partition uses n element breaks instead of 2n nodal intervals.
This is an isolated-process wrapper; physics and all Cartesian native columns remain.
"""
import geometry_compatible_ancf_stream as f
import stable_curvature_profiles as stable
import curvature_work_profiles as raw

_profiles=stable.profiles
_evaluate=raw.evaluate

def profiles(x,n,axis,kind,metric):
    assert n%2==0
    return _profiles(x,n//2,axis,kind,metric)

def evaluate(x,n,axis,metric):
    assert n%2==0
    return _evaluate(x,n//2,axis,metric)

stable.profiles=profiles
raw.evaluate=evaluate
