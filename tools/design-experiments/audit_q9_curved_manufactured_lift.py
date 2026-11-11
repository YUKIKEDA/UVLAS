"""Check the analytic (uninterpolated) manufactured field through thickness."""
from pathlib import Path
import json,itertools
import numpy as np
import check_q9_curved_bounded_lift as core
def values(name,a,b,x,y):
    r,d,c,n,_,_=core.geometry.surface(name,x,y);sx=np.sin(a*np.pi*(x+.5));sy=np.sin(b*np.pi*(y+.5));w=sx*sy;wx=a*np.pi*np.cos(a*np.pi*(x+.5))*sy;wy=b*np.pi*sx*np.cos(b*np.pi*(y+.5));_,_,_,_,rx,ry=core.geometry.surface(name,x,y);length=np.sqrt(d@d);dd=np.linalg.solve(np.array([rx,ry,d]),np.array([-length*wx,-length*wy,0*x]));k=.2 if name=='NORMAL-C' else 0.
    return r,d,c,w*n,dd,k*dd
def main():
    target=core.F/'q9-curved-manufactured-lift-audit.json'
    if target.exists():raise FileExistsError(target)
    cases=[]
    for name,(a,b) in itertools.product(['SADDLE','NORMAL-C'],core.WAVES):
        maximum=0.
        for x,y,z in itertools.product([-.31,-.09,.23],[-.27,.07,.34],[-.05,0.,.05]):
            r,d,c,u,v,t=values(name,a,b,x,y);dx=[v.imag/1e-25 for v in values(name,a,b,x+1e-25j,y)];dy=[v.imag/1e-25 for v in values(name,a,b,x,y+1e-25j)];gz=d+z*c;duz=v+z*t
            for derivative in [dx,dy]:
                ra,da,ca,ua,va,ta=derivative;ga=ra+z*da+z*z*ca/2;dua=ua+z*va+z*z*ta/2;maximum=max(maximum,float(abs(ga@duz+gz@dua)))
            maximum=max(maximum,float(abs(d@v)))
        cases.append(dict(geometry=name,wave=[a,b],maximum_covariant_shear_or_stretch=maximum,checks_pass=bool(maximum<=1e-11)))
    result=dict(source_sha256=core.digest(__file__),protocol_sha256=core.digest(core.F/'q9-curved-bounded-lift-experiment.md'),dependencies={Path(m.__file__).name:core.digest(m.__file__) for m in [core,core.geometry]},cases=cases,all_criteria_pass=bool(all(c['checks_pass'] for c in cases)),interpolated_field_certified=False)
    target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps(result))
    if not result['all_criteria_pass']:raise SystemExit(1)
if __name__=='__main__':main()
