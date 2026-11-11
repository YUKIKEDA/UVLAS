"""QR zero coefficients, general-LU complex derivative, tensor/H1 audit."""
from pathlib import Path
import hashlib,json,itertools,time
import numpy as np
from scipy.linalg import qr,solve,solve_triangular
import gsp_q3_shared_space as method
candidate=method.candidate
ROOT=Path(__file__).resolve().parents[2]
LAM=.3/(1.3*.4);MU=1/2.6
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rel(a,b):return float(np.linalg.norm(a)/max(1.,np.linalg.norm(b)))
def energy(E):return LAM*np.trace(E,axis1=1,axis2=2)**2+2*MU*np.einsum('pij,pij->p',E,E)

def coefficients(state,points):
    p,w=candidate.BASE.quadrature(9);mu,_=candidate.measure(state,0.,p,w);S=candidate.stress_field(state,p)
    _,E=candidate.linear_raw(state['q'],0.,p)
    A=(np.einsum('ij,pjn->pin',candidate.CIH,S)*np.sqrt(mu)[:,None,None]).reshape(-1,S.shape[2])
    Y=(np.einsum('ij,jk,pkn->pin',candidate.CIH,candidate.C,E)*np.sqrt(mu)[:,None,None]).reshape(-1,150)
    Q,R=qr(A,mode='economic',check_finite=False);a0=solve_triangular(R,Q.T@Y,check_finite=False)
    pc,wc,Vc,Sc,Hc=candidate.system(state,1e-25j)
    _,Ec=candidate.linear_raw(state['q'],1e-25j,pc)
    lc=np.einsum('p,pki,pkn->in',wc,Sc,Ec,optimize=True);ac=solve(Hc,lc,assume_a='gen',check_finite=False)
    field=candidate.stress_field(state,points)
    B0=np.einsum('ij,pjn->pin',candidate.CI,field@a0)[:,:,:144]
    B1=2*np.einsum('ij,pjn->pin',candidate.CI,field@(ac.imag/1e-25))[:,:,:144]
    _,J,_=candidate.kinematics(state['q'],0.,points)
    return B0,B1,J,a0,ac.imag/1e-25

def main():
    dest=ROOT/'docs/verification/gsp-q3-shared-audit.json';assert not dest.exists()
    start=time.monotonic();path=ROOT/'docs/verification/gsp-q3-shared-result.json';data=json.loads(path.read_text())
    for source,sha in data['sources'].items():assert digest(ROOT/source)==sha
    p,w=np.polynomial.legendre.leggauss(11);plane=np.array(list(itertools.product(p,p)));points=np.column_stack((plane,np.ones(len(plane))))
    weights=np.array([a*b for a,b in itertools.product(w,w)]);N,g=candidate.BASE.geometry.shape(plane)
    records=[];cache={};map_checks={};stop=None
    try:
        for record in data['records']:
            assert digest(ROOT/record['artifact'])==record['artifact_sha256'];saved=np.load(ROOT/record['artifact'])
            free=saved['free'];lookup={int(v):i for i,v in enumerate(free)}
            totals={kind:dict(error=0.,reference=0.,zero=0.,h1=0.) for kind in record['solutions']}
            kernel_square=matrix_square=area=control_zero=control_error=control_reference=0.
            for ids,key in zip(saved['ids'],record['map_keys']):
                if key not in cache:
                    item=data['maps'][key];assert digest(ROOT/item['artifact'])==item['artifact_sha256']
                    m=np.load(ROOT/item['artifact']);state={name:m[name] for name in method.KEYS}
                    B0,B1,J,a0,ad=coefficients(state,points);cache[key]=(B0,B1,J)
                    map_checks[key]=dict(zero=rel(a0-m['a0'],m['a0']),first=rel(ad-m['ad'],m['ad']))
                B0,B1,J=cache[key];weight=weights*np.linalg.det(J);K=np.zeros((144,saved['kernel'].shape[1]))
                for local,node in enumerate(ids):
                    if int(node) in lookup:K[9*local:9*local+6]=saved['kernel'][6*lookup[int(node)]:6*lookup[int(node)]+6]
                WK=np.einsum('ij,pjn->pin',candidate.CH,B0@K)
                WB=np.einsum('ij,pjn->pin',candidate.CH,B0[:,:,method.common.RD])
                kernel_square+=float(np.sum(weight[:,None,None]*WK**2));matrix_square+=float(np.sum(weight[:,None,None]*WB**2))
                x=(N@saved['coords'][ids,0])[:,0];k=record['kappa'];avec=np.column_stack((np.ones_like(x),np.zeros_like(x),k*x));D=1+k*k*x*x
                Jsurf=np.einsum('na,pnb->pab',saved['coords'][ids,0],g);G=Jsurf.transpose(0,2,1)@Jsurf
                aw=weights*np.sqrt(np.linalg.det(G));area+=aw.sum()
                target0=-avec[:,:,None]*avec[:,None,:]/D[:,None,None]
                control=saved['control_u'][ids].reshape(144)
                delta=np.einsum('kij,pk->pij',candidate.MB,B1@control)-target0
                control_error+=weight@energy(delta);control_reference+=weight@energy(target0)
                control_zero+=float(np.sum(weight[:,None]*(np.einsum('ij,pj->pi',candidate.CH,B0@control))**2))
                for kind in totals:
                    u=saved['u_'+kind][ids];target=target0*(np.ones_like(x) if kind=='quadratic' else x*x)[:,None,None]
                    E=np.einsum('kij,pk->pij',candidate.MB,B1@u.reshape(144))
                    totals[kind]['error']+=weight@energy(E-target);totals[kind]['reference']+=weight@energy(target)
                    totals[kind]['zero']+=float(np.sum(weight[:,None]*np.einsum('ij,pj->pi',candidate.CH,B0@u.reshape(144))**2))
                    val=N@u[:,3:9];du=np.einsum('na,pnb->pab',u[:,3:9],g)
                    totals[kind]['h1']+=aw@(np.sum(val*val,axis=1)+np.einsum('pai,pij,paj->p',du,np.linalg.inv(G),du))
            kernel=float(np.sqrt(kernel_square/matrix_square)/np.linalg.norm(saved['kernel']))
            c0=float(np.sqrt(control_zero)/max(1.,np.linalg.norm(saved['control_u'][free,:6])))
            c1=float(np.sqrt(control_error/control_reference))
            checks=dict(kernel=kernel<=1e-9,control_reproduction=max(abs(c0-record['control0']),abs(c1-record['control1']))<=1e-8,
                coefficients=all(max(map_checks[name].values())<=1e-9 for name in record['map_keys']))
            observations={}
            for kind,value in totals.items():
                error=float(np.sqrt(value['error']/value['reference']));h1=float(np.sqrt(value['h1']/area))
                zero=float(np.sqrt(value['zero']/matrix_square)/np.linalg.norm(saved['u_'+kind][free,:6]))
                difference=abs(error-record['solutions'][kind]['observed_error']);hdiff=abs(h1/record['solutions'][kind]['observed_h1']-1)
                checks[kind]=difference<=1e-8 and hdiff<=1e-8 and zero<=1e-9
                observations[kind]=dict(error=error,error_difference=difference,h1=h1,h1_difference=hdiff,zero_residual=zero)
            records.append(dict(kappa=k,amplitude=record['amplitude'],n=2,kernel_residual=kernel,control0=c0,control1=c1,
                observations=observations,primary_numerical_checks_pass=all(record['checks'].values()),checks={k:bool(v) for k,v in checks.items()}))
            print(records[-1],flush=True)
            if not all(checks.values()):stop='required_shared_audit_gate_failed';break
    except Exception as exc:stop=repr(exc)
    sources=[Path(__file__),Path(method.__file__),Path(candidate.__file__),path,ROOT/'docs/verification/gsp-q3-shared-protocol.md']
    dest.write_text(json.dumps(dict(records=records,map_checks=map_checks,
        all_reproduction_checks=bool(records) and len(records)==len(data['records']) and stop is None,
        primary_all_configurations_complete=data['all_checks_pass'],independent_global_resolve=False,physical_adoption=False,stop_reason=stop,
        sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start),indent=2),encoding='utf-8')
if __name__=='__main__':main()
