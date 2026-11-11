"""All-coordinate matrix audit of ME-1 C0 recovery, independent quadrature."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,time
import numpy as np
import moment_position_recovery as r
import check_moment_enriched_local_v2 as local
f=r.f;ROOT=rd_paths.ROOT
def norm(a):return float(np.linalg.norm(a))
def relative(a,b):return norm(a)/max(1.,norm(b))
def run(name,n,qs,h):
    start=time.monotonic();states=[f.prepare(q,h) for q in qs];R=r.Recovery(states,n);records=[]
    # Assemble compatible reference positions from the native coordinates.
    X=np.zeros((R.size,3));known=np.zeros(R.size,bool)
    for cell,s in enumerate(states):
        U=R.native[cell]
        # Native internal blocks are exact identity columns.
        ids=np.argmax(U[:s['ni']],axis=1);X[ids]=s['X0'][:s['ni']];known[ids]=True
    for item in R.facet.values():
        s=states[item['cell']];face=item['face'];o=item['offset'];X[o:o+face['n']]=s['X0'][face['offset']:face['offset']+face['n']];known[o:o+face['n']]=True
    assert known.all()
    rng=np.random.default_rng(20261104);Q=rng.normal(size=(R.size,3));v=rng.normal(size=(R.size,3));F=np.array([[1.01,.02,0],[0,.99,.01],[.01,0,1.02]])
    for cell,s in enumerate(states):
        p,w=f.volume_points(15,9);mu=w/8;A=R.evaluate(cell,p);V=r.modal(p,f.VID);U=R.native[cell];project=V.T@(mu[:,None]*A)
        moment=relative(project-U[:f.NV],U[:f.NV]);faceerrors=[];boundarywork=np.zeros((R.size,3))
        P=np.array([[1.2,.3,-.2],[.4,.7,.1],[-.1,.2,1.1]])
        for face in s['faces']:
            fp,fw=f.old.face_points(face['axis'],face['sign'],15,9);Af=R.evaluate(cell,fp);Mf=r.facevalue(face,fp);target=U[face['offset']:face['offset']+face['n']]
            faceerrors.append(relative(Mf.T@((fw/4)[:,None]*Af)-target,target))
            _,_,_,_,Cf=f.geometry(s['q'],h,fp,True);traction=(face['sign']*Cf[:,face['axis'],:])@P.T;boundarywork+=Af.T@(fw[:,None]*traction)
        physical,_,_,_,_=f.geometry(s['q'],h,p,True);affine=relative(A@(X@F.T)-physical@F.T,physical)
        constant=relative(A@np.ones(R.size)-np.ones(len(p)),np.ones(len(p))) if False else 0.
        # Native projection of a constant is not an all-ones coefficient vector.
        translation=np.zeros((R.size,3))
        for c,st in enumerate(states):
            ut=f.reduce(st,lambda p:np.tile([.2,-.1,.3],(len(p),1)));NN=R.native[c];ids=np.argmax(NN[:st['ni']],axis=1);translation[ids]=ut[:st['ni']]
        for item in R.facet.values():
            ut=f.reduce(states[item['cell']],lambda p:np.tile([.2,-.1,.3],(len(p),1)));fa=item['face'];o=item['offset'];translation[o:o+fa['n']]=ut[fa['offset']:fa['offset']+fa['n']]
        rigid=relative(A@translation-np.tile([.2,-.1,.3],(len(p),1)),A@translation)
        nativework=U.T@(f.boundary_map(s)@P.T);work=relative(boundarywork-nativework,nativework)
        # Polynomial evaluation derivatives by complex step; force transpose is
        # tested against independently differentiated scalar point-load work.
        samples=p[::max(1,len(p)//19)];As=R.evaluate(cell,samples);load=rng.normal(size=(len(samples),3));generalized=As.T@load
        eps=1e-25;dwork=np.sum((As@(Q+1j*eps*v)).imag/eps*load);virtual=abs(dwork-np.sum(v*generalized))/max(1.,abs(dwork))
        G=R.gradient(cell,samples);normal=np.cross((G@X)[:,0,:],(G@X)[:,1,:]);normal/=np.linalg.norm(normal,axis=1)[:,None]
        extra_ids=np.argmax(U[[75]],axis=1);extra=norm(A[:,extra_ids]);records.append(dict(cell=cell,volume_moments=moment,face_moments=max(faceerrors),affine=affine,translation=rigid,constant_P_work=work,virtual_work=float(virtual),extra_response=extra,normal_unit=float(np.max(abs(np.linalg.norm(normal,axis=1)-1)))))
    shared=[]
    t=np.linspace(-1,1,11);pp=np.array([[a,b] for a in t for b in t])
    for cell,s in enumerate(states):
        j,i=divmod(cell,n)
        for axis,neighbor in [(0,cell+1),(1,cell+n)]:
            if (axis==0 and i==n-1) or (axis==1 and j==n-1):continue
            a=R.evaluate(cell,r.points(axis,1,pp));b=R.evaluate(neighbor,r.points(axis,-1,pp));shared.append(relative(a-b,a))
    # Rebuild with distinct integration; compare physical map, not coefficients.
    RI=r.Recovery(states,n,15,9);independent=max(relative(R.evaluate(c,states[c]['p'])-RI.evaluate(c,states[c]['p']),R.evaluate(c,states[c]['p'])) for c in range(n*n))
    # Homogeneous clamped displacement uses no internal/TB data in its skeleton.
    RC=r.Recovery(states,n,clamp=True);dq=rng.normal(size=(RC.size,3))
    for key,item in RC.facet.items():
        if key[0]==0 and key[1]==0:dq[item['offset']:item['offset']+item['face']['n']]=0
    support=max(norm(RC.evaluate(c,r.points(0,-1,pp))@dq) for c in range(n*n) if c%n==0)
    fields=['volume_moments','face_moments','affine','translation','constant_P_work','virtual_work','normal_unit']
    passed=max([x[k] for x in records for k in fields]+shared+[support])<=1e-9 and independent<=1e-8 and min(x['extra_response'] for x in records)>1e-6
    return dict(name=name,n=n,h=h,global_scalar=R.size,records=records,shared=max(shared,default=0.),independent=independent,support=support,seconds=time.monotonic()-start,passed=bool(passed))
def main():
    out=ROOT/'rd/structure/main-shell/moment-enriched/records/moment-position-recovery-result.json';assert not out.exists();records=[];stops=[]
    for name,k,n,amp in [('flat',0.,2,0.),('curved-distorted',.4,2,.2)]:
        grid=f.m.mesh(k,n,amp);qs=[grid['ref'][ids].reshape(-1) for ids in grid['cells']]
        try:rec=run(name,n,qs,.001);records.append(rec);print(name,rec,flush=True)
        except Exception as exc:stops.append(dict(name=name,reason=repr(exc)));print('STOP',name,repr(exc),flush=True)
    try:rec=run('quartic-reference',1,[local.prior.generic()],.001);records.append(rec);print('quartic',rec,flush=True)
    except Exception as exc:stops.append(dict(name='quartic-reference',reason=repr(exc)));print('STOP quartic',repr(exc),flush=True)
    paths=[Path(__file__),Path(r.__file__),Path(f.__file__),ROOT/'rd/structure/main-shell/moment-enriched/records/moment-enriched-protocol.md']
    out.write_text(json.dumps(dict(records=records,stops=stops,passed=len(records)==3 and not stops and all(x['passed'] for x in records),physical_adoption=False,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}),indent=2),encoding='utf-8')
if __name__=='__main__':main()
