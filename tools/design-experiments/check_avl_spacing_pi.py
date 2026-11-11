"""Independent high-precision spacing validation before corrected symmetry runs."""
from pathlib import Path
from decimal import Decimal,localcontext
import hashlib,json
import numpy as np

ROOT=Path(__file__).resolve().parents[2];BASE=ROOT/'build/reference-avl';FOLDER=ROOT/'docs/verification'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def arctan_inverse(n,precision):
    with localcontext() as ctx:
        ctx.prec=precision+12;x=Decimal(1)/n;power=x;total=x;k=1
        while True:
            power*=-x*x;term=power/(2*k+1);total+=term
            if abs(term)<Decimal(10)**(-precision-8):return +total
            k+=1
def pi(precision):
    with localcontext() as ctx:
        ctx.prec=precision
        return +(16*arctan_inverse(5,precision)-4*arctan_inverse(239,precision))
def cosine(x,precision=80):
    with localcontext() as ctx:
        ctx.prec=precision+12;term=Decimal(1);value=term;k=1
        while True:
            term*=-x*x/((2*k-1)*(2*k));value+=term
            if abs(term)<Decimal(10)**(-precision-8):return +value
            k+=1
def main():
    low=pi(80);high=pi(160)
    with localcontext() as ctx:ctx.prec=160;pi_difference=abs(low-high)
    old=np.loadtxt(BASE/'spacing-old.txt',skiprows=1);new=np.loadtxt(BASE/'spacing-pi64.txt',skiprows=1)
    assert old.shape==new.shape==(132,3) and np.array_equal(old[:,:2],new[:,:2])
    records=[]
    with localcontext() as ctx:
        ctx.prec=80
        for n in [3,5,9,17,33,65]:
            a=old[old[:,0]==n,2];b=new[new[:,0]==n,2];errors=[]
            for k,x in enumerate(b):
                theta=high*Decimal(k)/Decimal(n-1);reference=(1-cosine(theta))/2
                errors.append(float(abs(Decimal.from_float(float(x))-reference)))
            records.append(dict(n=n,old_mirror_error=float(max(abs(a+a[::-1]-1))),new_mirror_error=float(max(abs(b+b[::-1]-1))),
                new_reference_error=max(errors),old_midpoint=float(a[n//2]),new_midpoint=float(b[n//2])))
    failed=json.loads((FOLDER/'avl-symmetry-equivalence-initial-result.json').read_text())
    observed=max(r['geometry_difference'] for r in failed['records'])
    predicted=.8*abs(records[0]['old_midpoint']-.5)
    envelope=.4*abs(float(Decimal.from_float(3.1415926535)-high))/2+1e-14
    checks=dict(pi_self=pi_difference<=Decimal('1e-75'),new_spacing=all(r['new_reference_error']<=5e-16 for r in records),
        new_symmetry=all(r['new_mirror_error']<=2e-15 for r in records),old_geometry_explained=abs(observed-predicted)<=1e-15)
    result=dict(protocol='docs/verification/avl-symmetry-spacing-pi-addendum.md',source_sha256=digest(Path(__file__)),
        pi80=str(low),pi160=str(high),pi_difference=str(pi_difference),records=records,observed_old_geometry_difference=observed,
        predicted_old_geometry_difference=predicted,old_new_coordinate_envelope=envelope,checks=checks,all_criteria_pass=all(checks.values()),
        files={p.name:digest(p) for p in [BASE/'spacing-old.txt',BASE/'spacing-pi64.txt']})
    (FOLDER/'avl-spacing-pi-result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(checks=checks,predicted=predicted,observed=observed,envelope=envelope),indent=2))

if __name__=='__main__':main()
