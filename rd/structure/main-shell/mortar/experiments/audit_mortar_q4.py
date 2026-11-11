"""Independent quadrature/basis audit; no independent global optimization."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,time
import numpy as np
import mortar_q4 as f
m=f.m;ROOT=rd_paths.ROOT
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def relative(a,b):return float(np.linalg.norm(a)/max(1e-30,np.linalg.norm(b)))
def energy(E):return .3/(1.3*.4)*np.trace(E,axis1=1,axis2=2)**2+2/2.6*np.einsum('pij,pij->p',E,E)

def main():
    source=ROOT/'rd/structure/main-shell/mortar/records/mortar-q4-result.json';parent=json.loads(source.read_text());dest=ROOT/'rd/structure/main-shell/mortar/records/mortar-q4-audit.json';assert not dest.exists()
    for p,h in parent['sources'].items():assert digest(ROOT/p)==h
    start=time.monotonic();locals=[];interfaces=[];records=[];stop=None
    try:
        for r in parent['local_records']:
            assert digest(ROOT/r['artifact'])==r['artifact_sha256'];s=np.load(ROOT/r['artifact'])
            U,force,tangent=f.raw_nonlinear(s['q'],s['current'],r['h'],s['v1'],s['v2'],11,7,True)
            obs=dict(energy=relative(U-s['energy'],s['energy']),force=relative(force-s['force'],s['force']),tangent=relative(tangent-s['tangent'],s['tangent']))
            locals.append(dict(geometry=r['geometry'],h=r['h'],observations=obs,passed=all(v<=1e-8 for v in obs.values())))
            if not locals[-1]['passed']:raise ValueError('local_independent_audit_failed')
        for r in parent['interfaces']:
            mesh=f.broken_mesh(r['kappa'],2,r['amplitude']);scans=[f.rigid_scan(mesh,h,r['degree'],True) for h in [.01,.0001]]
            same=all((a['ranks'],a['supported_ranks'])==(b['ranks'],b['supported_ranks']) for a,b in zip(scans,r['scans']))
            patches=[f.patch(mesh,.01,r['degree'],np.eye(3)+.03*b,True) for b in m.MB]
            checks=dict(rigid_rank_reproduced=same,patch=all(p['residual']<=1e-9 and p['gap']<=1e-9 for p in patches))
            interfaces.append(dict(kappa=r['kappa'],amplitude=r['amplitude'],degree=r['degree'],scans=scans,patches=patches,checks=checks))
            if not all(checks.values()):raise ValueError('interface_independent_audit_failed')
        for r in parent['records']:
            assert digest(ROOT/r['artifact'])==r['artifact_sha256'];saved=np.load(ROOT/r['artifact']);mesh=f.broken_mesh(r['kappa'],r['n'],r['amplitude'])
            p,w=m.quadrature(11,2);lookup={int(v):i for i,v in enumerate(saved['free_rd'])};totals={kind:dict(error=0.,reference=0.,h1=0.) for kind in r['solutions']};area=0.;kernel2=matrix2=0.
            for ids in mesh['cells']:
                q=mesh['ref'][ids].reshape(225);_,J,B0=m.raw(q,p,0.,True);_,_,Bc=m.raw(q,p,1e-25j,True);B1=Bc.imag/1e-25;weight=w*np.linalg.det(J)
                N=np.zeros((225,saved['kernel'].shape[1]))
                for a,node in enumerate(ids):
                    if int(node) in lookup:N[3*a:3*a+3]=saved['kernel'][3*lookup[int(node)]:3*lookup[int(node)]+3]
                kernel2+=np.sum(weight[:,None,None]*(B0@N)**2);matrix2+=np.sum(weight[:,None,None]*B0**2)
                Nr,gr=m.shape(p,0,True);J2=np.einsum('na,pnb->pab',mesh['ref'][ids[:25]],gr);metric=J2.transpose(0,2,1)@J2;aw=w*np.sqrt(np.linalg.det(metric));area+=aw.sum();x=(Nr@mesh['ref'][ids[:25]])[:,0]
                for kind,t in totals.items():
                    u=saved['u_'+kind][ids];observed=np.einsum('kij,pk->pij',m.MB,B1@u.reshape(-1));target=np.einsum('kij,pk->pij',m.MB,m.target(kind,r['kappa'],x))
                    t['error']+=weight@energy(observed-target);t['reference']+=weight@energy(target)
                    for field in [1,2]:
                        ns,g=m.shape(p,field,True);v=ns@u[m.SLICES[field]];du=np.einsum('na,pnb->pab',u[m.SLICES[field]],g)
                        t['h1']+=aw@(np.sum(v*v,axis=1)+np.einsum('pai,pij,paj->p',du,np.linalg.inv(metric),du))
            G,_=f.constraints(mesh,0.,r['degree'],11,7,True);Kr=np.kron(G[:,saved['free_rd']],np.eye(3));Kc=np.kron(G[:,saved['free_c']],np.eye(3))
            checks=dict(kernel=float(np.sqrt(kernel2/matrix2)/np.linalg.norm(saved['kernel']))<=1e-9,
                rd_constraint=relative(Kr@saved['kernel'],G)/max(1.,np.linalg.norm(saved['kernel']))<=1e-10,
                c_constraint=relative(Kc@saved['c_basis'],G)/max(1.,np.linalg.norm(saved['c_basis']))<=1e-10);observations={}
            for kind,t in totals.items():
                e=float(np.sqrt(t['error']/t['reference']));h1=float(np.sqrt(t['h1']/area));de=abs(e-r['solutions'][kind]['error']);dh=abs(h1/r['solutions'][kind]['h1']-1)
                observations[kind]=dict(error=e,error_difference=de,h1_difference=dh);checks[kind]=de<=1e-8 and dh<=1e-8
            records.append(dict(kappa=r['kappa'],amplitude=r['amplitude'],degree=r['degree'],n=r['n'],observations=observations,checks={k:bool(v) for k,v in checks.items()}))
            if not all(checks.values()):raise ValueError('shared_independent_audit_failed')
    except Exception as exc:stop=repr(exc)
    sources=[Path(__file__),Path(f.__file__),Path(m.__file__),source,ROOT/'rd/structure/main-shell/mortar/records/mortar-q4-protocol.md']
    dest.write_text(json.dumps(dict(local_records=locals,interfaces=interfaces,records=records,
        all_reproduction_checks=len(locals)==len(parent['local_records']) and len(interfaces)==len(parent['interfaces']) and len(records)==len(parent['records']) and stop is None,
        stop_reason=stop,independent_global_resolve=False,physical_adoption=False,sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start),indent=2),encoding='utf-8')
    print('audit',len(locals),len(interfaces),len(records),'stop',stop)
if __name__=='__main__':main()
