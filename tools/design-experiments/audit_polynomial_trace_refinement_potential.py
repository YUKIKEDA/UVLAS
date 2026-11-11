"""Physical tensor/jump-energy observation of stored globally coupled trace fits."""
from pathlib import Path
import hashlib,json,time
import numpy as np
import polynomial_trace_ancf as f
m=f.m;ROOT=Path(__file__).resolve().parents[2]
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def energy(E):return .3/(1.3*.4)*np.trace(E,axis1=1,axis2=2)**2+2/2.6*np.einsum('pij,pij->p',E,E)
def main():
    source=ROOT/'docs/verification/polynomial-trace-potential-refinement-result.json';parent=json.loads(source.read_text());dest=ROOT/'docs/verification/polynomial-trace-potential-refinement-audit.json';assert not dest.exists()
    start=time.monotonic();records=[];stops=[]
    for r in parent['records']:
        try:
            assert digest(ROOT/r['artifact'])==r['artifact_sha256'];s=np.load(ROOT/r['artifact']);err=np.zeros(2);ref=np.zeros(2);H=np.zeros(2);area=0.
            for c,ids in enumerate(s['cells']):
                q=s['ref'][ids].reshape(-1);values=s['local'][c]
                if r['kind']=='baseline':
                    p,w=f.volume_points(15,9);_,g,_=m.tables(p,r['h']*p[:,2]/2,True);J=np.einsum('na,pnb->pab',q.reshape(75,3),g);G=g@np.linalg.inv(J);weight=w*np.linalg.det(J)*r['h']/2;state=None
                else:state=f.prepare(q,r['h'],r['kind'],15,9,True);p=state['p'];weight=state['weight'];G=state['Gw'].transpose(0,2,1)
                N,_=m.shape(p,0,True);x=(N@q.reshape(75,3)[:25])[:,0];z=r['h']*p[:,2]/2;a=np.column_stack((np.ones(len(x)),np.zeros(len(x)),r['kappa']*x));den=1+r['kappa']**2*x*x-r['kappa']*z
                plane,wp=m.quadrature(15,2);Ns,gs=m.shape(plane,0,True);Js=np.einsum('na,pnb->pab',q.reshape(75,3)[:25],gs);metric=Js.transpose(0,2,1)@Js;aw=wp*np.sqrt(np.linalg.det(metric));area+=aw.sum()
                for j,kind in enumerate(['quartic','exponential']):
                    u=values[:,j].reshape(-1,3);du=np.einsum('na,pnj->paj',u,G);observed=(du+du.transpose(0,2,1))/2
                    curvature=x*x if kind=='quartic' else np.exp(x);target=-z[:,None,None]*curvature[:,None,None]/den[:,None,None]*a[:,:,None]*a[:,None,:]
                    err[j]+=weight@energy(observed-target);ref[j]+=weight@energy(target)
                    if state is not None:
                        for face in state['faces']:err[j]+=state['tau']*face['measure']@np.sum((face['jump']@u)**2,axis=1)
                    for field in [1,2]:
                        nf,gf=m.shape(plane,field,True);uv=u[m.SLICES[field]];v=nf@uv;dv=np.einsum('na,pnj->paj',uv,gf)
                        H[j]+=aw@(np.sum(v*v,axis=1)+np.einsum('pai,pij,paj->p',dv,np.linalg.inv(metric),dv))
            observed_error=np.sqrt(err/ref);observed_H=np.sqrt(H/area);de=max(abs(x-y) for x,y in zip(observed_error,r['observed']['error']));dh=max(abs(x/y-1) for x,y in zip(observed_H,r['observed']['h1']))
            records.append(dict(kind=r['kind'],kappa=r['kappa'],amplitude=r['amplitude'],h=r['h'],error=observed_error.tolist(),h1=observed_H.tolist(),error_difference=de,h1_difference=dh,passed=bool(de<=1e-8 and dh<=1e-8)))
        except Exception as exc:stops.append(dict(kind=r['kind'],kappa=r['kappa'],amplitude=r['amplitude'],h=r['h'],reason=repr(exc)))
    sources=[Path(__file__),Path(f.__file__),source,ROOT/'docs/verification/polynomial-trace-protocol.md']
    dest.write_text(json.dumps(dict(records=records,stops=stops,all_reproduction_checks=len(records)==len(parent['records']) and all(r['passed'] for r in records) and not stops,
        independent_global_resolve=False,physical_adoption=False,sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start),indent=2),encoding='utf-8')
    print('shared audit',len(records),'failed',[r for r in records if not r['passed']],'stops',stops)
if __name__=='__main__':main()
