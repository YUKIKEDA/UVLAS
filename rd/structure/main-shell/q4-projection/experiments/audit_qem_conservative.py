import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,itertools,time
import numpy as np
from scipy.linalg import svd
import qem_el2_q4 as m
ROOT=rd_paths.ROOT
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rel(a,b):return float(np.linalg.norm(a)/max(1.,np.linalg.norm(b)))
def energy(E):return .3/(1.3*.4)*np.trace(E,axis1=1,axis2=2)**2+2/2.6*np.einsum('pij,pij->p',E,E)
def fields(q,p,method):
    _,J,B0=m.raw(q,p,0.,True);_,_,Bc=m.raw(q,p,1e-25j,True);B1=Bc.imag/1e-25
    if method=='raw':return J,B0,B1
    pn,_=m.native_plane();_,_,Bn=m.raw(q,pn,0.,True);_,_,Bcn=m.raw(q,pn,1e-25j,True);N,_=m.shape(p,0,True)
    R0=np.einsum('pn,nki->pki',N,Bn);R1=np.einsum('pn,nki->pki',N,Bcn.imag/1e-25)
    if method=='conservative':
        pg,wg=m.quadrature(11,2);_,Jg,Bg=m.raw(q,pg,0.,True);Ng,_=m.shape(pg,0,True);mu=wg*np.linalg.det(Jg);mu/=mu.sum()
        R0+=np.einsum('p,pki->ki',mu,Bg-np.einsum('pn,nki->pki',Ng,Bn))
    return J,R0,R1

def independent_force(q,current,h):
    plane,wp=m.quadrature(11,2);z,wz=np.polynomial.legendre.leggauss(7)
    points=np.array([[s,t,c] for (s,t),c in itertools.product(plane,z)]);w=np.array([a*b for a,b in itertools.product(wp,wz)])
    pn,_=m.native_plane();native=np.array([[s,t,c] for (s,t),c in itertools.product(pn,z)])
    def data(p):
        _,g,_=m.tables(p,h*p[:,2]/2,True);J=np.einsum('na,pnb->pab',q.reshape(75,3),g);G=g@np.linalg.inv(J)
        F=np.eye(3)+np.einsum('na,pnb->pab',(current-q).reshape(75,3),G)
        E=np.einsum('kij,pij->pk',m.MB,(F.transpose(0,2,1)@F-np.eye(3))/2)
        B=np.einsum('kij,pai,pnj->pkna',m.MB,F,G).reshape(len(p),6,225)
        return J,E,B
    J,E,B=data(points);Jn,En,Bn=data(native);mu=w*np.linalg.det(J);mu/=mu.sum();N,_=m.shape(plane,0,True)
    IE=np.einsum('pn,nzk->pzk',N,En.reshape(25,7,6)).reshape(-1,6)
    Er=IE+np.einsum('p,pk->k',mu,E-IE);stress=Er@m.C.T;mean=np.einsum('p,pk->k',mu,stress)
    coef=np.einsum('pz,pn,pzk->nzk',mu.reshape(121,7),N,(stress-mean).reshape(121,7,6))
    force=np.einsum('p,pki,k->i',mu,B,mean)+np.einsum('nzki,nzk->i',Bn.reshape(25,7,6,225),coef)
    W=np.einsum('p,pk,pk->',mu,Er,stress)/2
    return W,force

def main():
    path=ROOT/'rd/structure/main-shell/q4-projection/records/qem-conservative-result.json';parent=json.loads(path.read_text());dest=ROOT/'rd/structure/main-shell/q4-projection/records/qem-conservative-audit.json';assert not dest.exists()
    for p,h in parent['sources'].items():assert digest(ROOT/p)==h
    start=time.monotonic();locals=[];records=[];stop=None;p,w=m.quadrature(11,2);Nr,gr=m.shape(p,0,True)
    try:
        for r in parent['local_records']:
            assert digest(ROOT/r['artifact'])==r['artifact_sha256'];s=np.load(ROOT/r['artifact'])
            U,f=independent_force(s['q'],s['current'],r['h']);Uc,fc=independent_force(s['q'],s['current']+1e-25j*s['v2'],r['h'])
            metrics=dict(energy=rel(U-s['energy'],s['energy']),force=rel(f-s['force'],s['force']),
                tangent=rel(s['v1']@(fc.imag/1e-25)-s['tangent'],s['tangent']))
            locals.append(dict(geometry=r['geometry'],h=r['h'],metrics=metrics,passed=all(v<=1e-8 for v in metrics.values())))
            if not locals[-1]['passed']:raise RuntimeError('independent_dual_force_audit_failed')
        for r in parent['records']:
            assert digest(ROOT/r['artifact'])==r['artifact_sha256'];s=np.load(ROOT/r['artifact']);free=s['free_rd'];lookup={int(v):i for i,v in enumerate(free)}
            totals={kind:dict(error=0.,reference=0.,h1=0.,zero=0.) for kind in r['solutions']};area=kernel2=matrix2=0.
            for ids in s['cells']:
                q=s['ref'][ids].reshape(225);J,B0,B1=fields(q,p,r['method']);weight=w*np.linalg.det(J);K=np.zeros((225,s['kernel'].shape[1]))
                for local,node in enumerate(ids[:50]):
                    if int(node) in lookup:K[3*local:3*local+3]=s['kernel'][3*lookup[int(node)]:3*lookup[int(node)]+3]
                WB=np.einsum('ij,pjn->pin',m.CH,B0[:,:,:150]);WK=np.einsum('ij,pjn->pin',m.CH,B0@K)
                matrix2+=np.sum(weight[:,None,None]*WB**2);kernel2+=np.sum(weight[:,None,None]*WK**2)
                surface=np.einsum('na,pnb->pab',s['ref'][ids[:25]],gr);G=surface.transpose(0,2,1)@surface;aw=w*np.sqrt(np.linalg.det(G));area+=aw.sum();x=(Nr@s['ref'][ids[:25]])[:,0]
                k=r['kappa'];a=np.column_stack((np.ones_like(x),np.zeros_like(x),k*x));D=1+k*k*x*x
                for kind in totals:
                    v={'quadratic':np.ones_like(x),'cubic':x,'quartic':x*x,'exponential':np.exp(x)}[kind];target=-v[:,None,None]*a[:,:,None]*a[:,None,:]/D[:,None,None]
                    u=s['u_'+kind][ids];E=np.einsum('kij,pk->pij',m.MB,B1@u.reshape(225))
                    totals[kind]['error']+=weight@energy(E-target);totals[kind]['reference']+=weight@energy(target)
                    totals[kind]['zero']+=np.sum(weight[:,None]*np.einsum('ij,pj->pi',m.CH,B0@u.reshape(225))**2)
                    for field in [1,2]:
                        N,g=m.shape(p,field,True);value=N@u[m.SLICES[field]];du=np.einsum('na,pnb->pab',u[m.SLICES[field]],g)
                        totals[kind]['h1']+=aw@(np.sum(value*value,axis=1)+np.einsum('pai,pij,paj->p',du,np.linalg.inv(G),du))
            kernel=float(np.sqrt(kernel2/matrix2)/np.linalg.norm(s['kernel']));checks=dict(kernel=kernel<=1e-9);observed={}
            for kind,t in totals.items():
                e=float(np.sqrt(t['error']/t['reference']));h1=float(np.sqrt(t['h1']/area))
                diff=abs(e-r['solutions'][kind]['error']);hdiff=abs(h1/r['solutions'][kind]['h1']-1)
                observed[kind]=dict(error=e,error_difference=diff,h1_difference=hdiff);checks[kind]=diff<=1e-8 and hdiff<=1e-8
            records.append(dict(method=r['method'],kappa=r['kappa'],amplitude=r['amplitude'],n=r['n'],observations=observed,checks={k:bool(v) for k,v in checks.items()}))
            if not all(checks.values()):raise RuntimeError('shared_field_audit_failed')
    except Exception as exc:stop=repr(exc)
    sources=[Path(__file__),Path(m.__file__),path,ROOT/'rd/structure/main-shell/q4-projection/records/qem-conservative-protocol.md']
    dest.write_text(json.dumps(dict(local_records=locals,records=records,
        all_reproduction_checks=bool(locals) and len(locals)==len(parent['local_records']) and len(records)==len(parent['records']) and stop is None,
        stop_reason=stop,independent_global_resolve=False,physical_adoption=False,
        sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start),indent=2),encoding='utf-8')
    print('local audits',len(locals),'shared audits',len(records),'stop',stop)
if __name__=='__main__':main()
