"""Coefficient identities and unrestricted local thin constraints, nine cases."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,time,itertools
import numpy as np
from scipy.linalg import svd
import reference_work_thin_limit as t
import check_reference_work_local as prior
import check_q9_plate_convergence as memory
ROOT=rd_paths.ROOT
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def relative(a,b):return float(np.linalg.norm(a-b)/max(1.,np.linalg.norm(b)))
def reference_series(op):
    s=op['state'];_,modal=t.c.geometry_coefficients(np.asarray(s['q'],dtype=np.float64).tobytes());out=[]
    for power in range(3):
        vol=np.zeros((75,3))
        for j,(i,k,z) in enumerate(t.f.VID):
            norm=np.sqrt((2*i+1)*(2*k+1)*(2*z+1))
            if power==0 and z==0:vol[j]=modal[0,i,k]/norm
            if power==1 and z==1:vol[j]=modal[1,i,k]/(2*norm)
            if power==2 and z==0:vol[j]=modal[2,i,k]/(24*norm)
            if power==2 and z==2:vol[j]=modal[2,i,k]/(12*norm)
        coeff=np.zeros((len(t.c.IDS),3));coeff[t.f.VI]=vol;value=np.zeros_like(s['X0']);value[:75]=vol
        for face in s['faces']:value[face['offset']:face['offset']+face['n']]=face['Mcoeff'].T@t.c.restrict(coeff,face['axis'],face['sign'])
        out.append(value)
    return out
def evaluate(q,h,name):
    op=t.operators(q);s=op['state'];finite=t.f.prepare(q,h);nc=3*s['nscalar'];C=sum(h**i*op['C'][:,:,:,i] for i in range(5));G=np.einsum('pjn,pjk->pnk',s['Gwref'],C,optimize=True)
    A=np.einsum('kab,pnb->pkna',t.m.MB,G).reshape(len(op['p']),6,nc);main=relative(A,finite['det'][:,None,None]*finite['B'])
    si=t.f.prepare(q,h,order=15,zorder=9,independent=True);pi=si['p'];L,_=t.c.scalar_basis(pi,True);cf,_=t.c.geometry_coefficients(np.asarray(q,dtype=np.float64).tobytes(),True);Ci=np.einsum('pm,mabj->pabj',L,cf.reshape(len(t.c.IDS),3,3,5),optimize=True)
    body=np.einsum('pm,mjr->pjr',L,s['bodycoeff'],optimize=True)@s['D'];Gh=np.einsum('pjn,pjk->pnk',body,sum(h**i*Ci[:,:,:,i] for i in range(5)),optimize=True)
    Ai=np.einsum('kab,pnb->pkna',t.m.MB,Gh).reshape(len(pi),6,nc);ind=relative(Ai,si['det'][:,None,None]*si['B'])
    J=np.vstack([x['hat'].T@((x['w']/4)[:,None]*x['jump']) for x in s['faces']]);JTB=J[:sum(x['n'] for x in s['faces'][:2])];Z0,info=t.nullspace(np.vstack((J,s['Gwref'][:,2,:])))
    ids=[v for v in t.f.PID if v[2]==0 and v[:2]!=(5,5)];coeff=np.zeros((len(t.c.IDS),35))
    for i,v in enumerate(ids):coeff[t.c.LOOK[v],i]=1.
    P=np.zeros((s['nscalar'],35));P[:75]=coeff[t.f.VI]
    for face in s['faces']:P[face['offset']:face['offset']+face['n']]=face['Mcoeff'].T@t.c.restrict(coeff,face['axis'],face['sign'])
    Q=np.linalg.qr(P,mode='reduced')[0];span=relative(Z0-Q@(Q.T@Z0),Z0)
    roots,_=t.leading_rows(op,0.);series=reference_series(op);ridx=[0]+[face['offset'] for face in s['faces']];rigid=[]
    for axis in np.eye(3):
        trans=np.zeros_like(s['X0']);trans[ridx]=axis;rigid.append([trans,np.zeros_like(trans),np.zeros_like(trans)]);rigid.append([np.cross(axis,x) for x in series])
    C01=t.constraints(op);cr=er=0.
    for q0,q1,q2 in rigid:
        pair=np.r_[q0.reshape(-1),q1.reshape(-1)];cr=max(cr,float(np.linalg.norm(C01@pair)/max(1.,np.linalg.norm(C01)*np.linalg.norm(pair))))
        E=roots[0]@q0.reshape(-1)+roots[1]@q1.reshape(-1)+roots[2]@q2.reshape(-1);er=max(er,float(np.linalg.norm(E)/max(1.,sum(np.linalg.norm(x) for x in roots))))
    ranks_J=t.rank_info(svd(J,compute_uv=False));ranks_TB=t.rank_info(svd(JTB,compute_uv=False))
    checks=dict(coefficient_identity=main<=1e-9,independent_identity=ind<=1e-8,stabilization_kernel=ranks_J==[s['nscalar']-105]*3,
        leading_kernel=info['nullity']==35 and span<=1e-9,rigid=cr<=1e-9 and er<=1e-9)
    return dict(geometry=name,h=h,coordinate_count=nc,coefficient_identity=main,independent_identity=ind,J_ranks=ranks_J,TB_ranks=ranks_TB,leading_kernel=info,polynomial_span_difference=span,
        rigid_constraint=cr,rigid_bending=er,checks={k:bool(v) for k,v in checks.items()})
def main():
    dest=ROOT/'rd/structure/main-shell/reference-work/records/reference-work-thin-local-result.json';assert not dest.exists();start=time.monotonic();records=[];stops=[]
    qs=[('flat',t.m.mesh(0.,1,0.)),('curved-distorted',t.m.mesh(.4,2,.2))];geoms=[(name,g['ref'][g['cells'][0]].reshape(-1)) for name,g in qs]+[('quartic-reference',prior.generic())]
    for name,q in geoms:
        for h in [.01,.001,.0001]:
            try:
                tick=time.monotonic();r=evaluate(q,h,name);r['checks']['resources']=time.monotonic()-tick<=900 and memory.memory_bytes()<=2*1024**3;records.append(r);print(name,h,r['checks'],flush=True)
            except Exception as exc:stops.append(dict(geometry=name,h=h,reason=repr(exc)));print('STOP',name,h,repr(exc),flush=True)
    paths=[Path(__file__),ROOT/'rd/structure/main-shell/reference-work/experiments/reference_work_thin_limit.py',ROOT/'rd/structure/main-shell/reference-work/records/reference-work-thin-protocol.md',ROOT/'rd/structure/main-shell/reference-work/records/reference-work-thin-preparation.json']
    dest.write_text(json.dumps(dict(records=records,stops=stops,all_checks_pass=len(records)==9 and not stops and all(all(r['checks'].values()) for r in records),main_shell_adopted=False,
        sources={p.relative_to(ROOT).as_posix():digest(p) for p in paths},seconds=time.monotonic()-start),indent=2),encoding='utf-8')
if __name__=='__main__':main()
