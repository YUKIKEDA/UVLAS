"""Exact interface rank of a specified pure-bending sector, tensor bound audit."""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import hashlib,json,itertools
import numpy as np
import mixed_order_150 as m
ROOT=Path(__file__).resolve().parents[2]
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rank(matrix):
    A=[[F(x) for x in row] for row in matrix];r=0
    for col in range(len(A[0])):
        pivot=next((i for i in range(r,len(A)) if A[i][col]),None)
        if pivot is None:continue
        A[r],A[pivot]=A[pivot],A[r];factor=A[r][col];A[r]=[v/factor for v in A[r]]
        for i in range(r+1,len(A)):
            if A[i][col]:
                factor=A[i][col];A[i]=[x-factor*y for x,y in zip(A[i],A[r])]
        r+=1
        if r==len(A):break
    return r
def interface_matrix(amp):
    mesh=m.mesh(.4,2,amp);k=F(2,5);rows=[];edges=[]
    # Unknowns per cell: integration constant Uc, W0,W1,W2,W3.
    fields=[{0:{0:F(1)},2:{2:-k/2},3:{3:-2*k/3},4:{4:-3*k/4}},
        {i+1:{i:F(1)} for i in range(4)},{i+1:{i-1:F(i)} for i in range(1,4)}]
    for left,right in itertools.combinations(range(4),2):
        shared=sorted(set(mesh['cells'][left,:25])&set(mesh['cells'][right,:25]))
        if len(shared)!=5:continue
        a=F.from_float(float(mesh['ref'][shared[0],0]));b=F.from_float(float(mesh['ref'][shared[-1],0]))-a
        edges.append(dict(cells=[left,right],a=str(a),b=str(b),x_varies=bool(b)))
        for field in fields:
            for power in range(5):
                row=[F(0)]*20
                for local,polynomial in field.items():
                    value=sum(coeff*comb(degree,power)*a**(degree-power)*b**power for degree,coeff in polynomial.items() if degree>=power)
                    row[5*left+local]+=value;row[5*right+local]-=value
                if any(row):rows.append(row)
    for cell,ids in enumerate(mesh['cells']):
        if np.count_nonzero(mesh['ref'][ids[:25],0]==0.)!=5:continue
        for local in [0,1,2]:
            row=[F(0)]*20;row[5*cell+local]=F(1);rows.append(row)
    return rows,edges
def main():
    path=ROOT/'docs/verification/mixed-150-root-diagnostic.json';parent=json.loads(path.read_text());dest=ROOT/'docs/verification/mixed-150-root-audit.json';assert not dest.exists()
    for file,sha in parent['sources'].items():assert digest(ROOT/file)==sha
    records=[]
    for amp,expected in [(0.,16),(.2,18)]:
        rows,edges=interface_matrix(amp);r=rank(rows)
        records.append(dict(amplitude=amp,rows=len(rows),unknowns=20,exact_rank=r,sector_dimension=20-r,expected_rank=expected,
            edges=edges,passed=r==expected,scope='ideal_bilinear_geometry_y_independent_cubic_W_quartic_U_sector'))
    bound=parent['exact_sector_bound'];alpha=float(F(bound['alpha']));beta=float(F(bound['beta']));k=.4;nu=.3
    p,w=np.polynomial.legendre.leggauss(64);x=(p+1)/2;w=w/2;D=1+k*k*x*x
    tangent=np.column_stack((np.ones_like(x),np.zeros_like(x),k*x))/np.sqrt(D)[:,None]
    normal=np.column_stack((-k*x,np.zeros_like(x),np.ones_like(x)))/np.sqrt(D)[:,None]
    tt=tangent[:,:,None]*tangent[:,None,:];nn=normal[:,:,None]*normal[:,None,:]
    delta=alpha+beta*x-x*x;res=-delta[:,None,None]*tt+(nu/(1-nu)*delta)[:,None,None]*nn;target=-x[:,None,None]**2*tt
    lam=nu/((1+nu)*(1-2*nu));mu=1/(2*(1+nu))
    energy=lambda E:lam*np.trace(E,axis1=1,axis2=2)**2+2*mu*np.einsum('pij,pij->p',E,E)
    observed=float((w*D)@energy(res)/((w*D)@energy(target)))
    exact=float(F(bound['relaxed_relative_error_squared']));difference=abs(observed-exact)
    stationarity=max(abs((w*D*delta)@np.ones_like(x)),abs((w*D*delta)@x))
    sources=[Path(__file__),Path(m.__file__),path,ROOT/'docs/verification/mixed-150-root-audit-protocol.md']
    dest.write_text(json.dumps(dict(interface_records=records,tensor_bound_squared=observed,exact_bound_squared=exact,bound_difference=difference,
        stationarity=float(stationarity),all_checks_pass=bool(all(r['passed'] for r in records) and max(difference,stationarity)<=1e-10),
        full_discrete_kernel_exactly_certified=False,finite_h_lower_bound=False,
        sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources}),indent=2),encoding='utf-8')
    print(json.dumps(dict(ranks=[r['exact_rank'] for r in records],sector_dimensions=[r['sector_dimension'] for r in records],bound_difference=difference,stationarity=stationarity),indent=2))
if __name__=='__main__':main()
