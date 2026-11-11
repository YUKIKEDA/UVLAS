"""Direct Decimal point formula and Clenshaw check of saved polynomial fields."""
from pathlib import Path
from decimal import Decimal,localcontext
import hashlib,json,itertools,time
import numpy as np
import gsp_q3_basis_v2 as v2
ROOT=Path(__file__).resolve().parents[2]
D=v2.D
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rel(a,b):return float(np.linalg.norm(a)/max(1.,np.linalg.norm(b)))

def direct(state,point):
    x,y,z=map(D,point);M=v2.old.BASE.geometry.COEFFICIENT
    def cardinal(t):
        return ([sum(D(M[p,i])*t**p for p in range(4)) for i in range(4)],
            [sum(D(M[p,i])*p*t**(p-1) for p in range(1,4)) for i in range(4)])
    ls,ds=cardinal(x);lt,dt=cardinal(y)
    N=[lt[j]*ls[i] for j,i in itertools.product(range(4),repeat=2)]
    gs=[lt[j]*ds[i] for j,i in itertools.product(range(4),repeat=2)]
    gt=[dt[j]*ls[i] for j,i in itertools.product(range(4),repeat=2)]
    q=state['q'].reshape(16,3,3)
    r=[[sum(N[n]*D(q[n,k,a]) for n in range(16)) for a in range(3)] for k in range(3)]
    A=[[sum(grad[n]*D(q[n,0,a]) for n in range(16))/D(state['frame_scale'][b]) if b<2 else r[1][a]/D(state['frame_scale'][2])
        for b,grad in enumerate([gs,gt,N])] for a in range(3)]
    det=sum(Decimal((-1)**sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3)))*A[0][perm[0]]*A[1][perm[1]]*A[2][perm[2]] for perm in itertools.permutations(range(3)))
    mapped=[]
    for k in range(6):
        tensor=[[sum(D(v2.old.MB[k,i,j])*A[a][i]*A[b][j] for i,j in itertools.product(range(3),repeat=2))/det for b in range(3)] for a in range(3)]
        mapped.append([sum(D(v2.old.MB[c,a,b])*tensor[a][b] for a,b in itertools.product(range(3),repeat=2)) for c in range(6)])
    legs=[]
    for t in [x,y,z]:legs.append([p*D(np.sqrt((2*k+1)/2)) for k,p in enumerate([Decimal(1),t,(3*t*t-1)/2,(5*t**3-3*t)/2])])
    columns=[]
    for c,imax,jmax,kmax in [(0,2,3,2),(1,3,2,2),(2,3,3,1),(3,3,2,2),(4,2,3,2),(5,2,2,2)]:
        for i,j,k in itertools.product(range(imax+1),range(jmax+1),range(kmax+1)):
            value=legs[0][i]*legs[1][j]*legs[2][k];columns.append([value*a for a in mapped[c]])
    phi=[Decimal(1)]+[(r[0][a]-D(state['center'].reshape(3)[a]))/D(state['length']) for a in range(3)]
    phi += [z*r[1][a] for a in range(3)]+[D(state['length'])*z*z*r[2][a] for a in range(3)]
    for p in phi:
        for c in range(6):columns.append([p if i==c else Decimal(0) for i in range(6)])
    T=np.array([[D(x) for x in row] for row in state['T']],dtype=object)
    return np.asarray(np.array(columns,dtype=object).T@T,dtype=float)

def main():
    dest=ROOT/'docs/verification/gsp-q3-basis-v2-audit.json';assert not dest.exists()
    start=time.monotonic();path=ROOT/'docs/verification/gsp-q3-basis-v2-result.json';parent=json.loads(path.read_text())
    for f,h in parent['sources'].items():assert digest(ROOT/f)==h
    records=[]
    for item in parent['records']:
        assert digest(ROOT/item['artifact'])==item['artifact_sha256']
        data=np.load(ROOT/item['artifact']);state={k:data[k] for k in ['q','frame_scale','length','center','T','active','column_norms']}
        points=data['arbitrary_points'];poly=dict(coeff=data['coeff120'],denominator=data['denominator'])
        clenshaw=np.polynomial.legendre.legval3d(points[:,0],points[:,1],points[:,2],poly['coeff']).transpose(2,0,1)
        den=np.polynomial.polynomial.polyval2d(points[:,0],points[:,1],poly['denominator'])
        clenshaw/=den[:,None,None]
        explicit=v2.evaluate(poly,points)
        with localcontext() as ctx:
            ctx.prec=120;expected=np.array([direct(state,p) for p in points])
        raw=v2.evaluate(dict(coeff=data['raw_coeff120'],denominator=data['denominator']),points)
        metrics=dict(direct_decimal_field=rel(raw-expected,expected),clenshaw_field=rel(explicit-clenshaw,explicit))
        grams=[]
        for order in [9,11,13]:
            p,w=v2.old.BASE.quadrature(order);mu,_=v2.old.measure(state,0.,p,w)
            S=np.polynomial.legendre.legval3d(p[:,0],p[:,1],p[:,2],poly['coeff']).transpose(2,0,1)
            S/=np.polynomial.polynomial.polyval2d(p[:,0],p[:,1],poly['denominator'])[:,None,None]
            H=np.einsum('p,pki,kl,plj->ij',mu,S,v2.old.CI,S,optimize=True)
            grams.append(dict(order=order,identity_error=rel(H-np.eye(len(H)),np.eye(len(H))),
                primary_difference=rel(H-data[f'gram{order}'],H)))
        checks=dict(fields=all(v<=1e-10 for v in metrics.values()),
            gram_reproduction=all(g['primary_difference']<=1e-10 for g in grams))
        records.append(dict(n=item['n'],cell=item['cell'],metrics=metrics,grams=grams,checks=checks))
        print(records[-1],flush=True)
    sources=[Path(__file__),Path(v2.__file__),path]
    dest.write_text(json.dumps(dict(records=records,all_reproduction_checks=all(all(r['checks'].values()) for r in records) and len(records)==len(parent['records']),
        primary_complete=parent['all_checks_pass'],physical_adoption=False,seconds=time.monotonic()-start,
        sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources}),indent=2),encoding='utf-8')
if __name__=='__main__':main()
