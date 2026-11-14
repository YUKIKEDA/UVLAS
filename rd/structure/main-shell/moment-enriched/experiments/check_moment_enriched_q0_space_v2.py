"""Post-derived full q0 image: test six-field hypothesis without imposing it."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,time
import numpy as np
from scipy.linalg import svd,qr
import check_moment_enriched_thin_global_v2 as solver
import moment_enriched_exact_isometry_v2 as exact
import check_reference_work_thin_q0_space_v2 as previous
import check_q9_plate_convergence as memory
t=solver.t;ROOT=rd_paths.ROOT
def field(q,p,k):
    old=previous.field(q,p,k);N,_=t.m.shape(p,0);x=(N@q.reshape(75,3)[:25])[:,0]
    return np.column_stack((old,-k*x**5/15,np.zeros(len(x)),x**4/12))
def span(a,b,rank):
    U,sv,_=svd(a,full_matrices=False,check_finite=False);ranks=t.rank_info(sv);Q=qr(b,mode='economic',check_finite=False)[0]
    qa=U[:,:ranks[1]];return dict(ranks=ranks,forward=float(np.linalg.norm(qa-Q@(Q.T@qa))),reverse=float(np.linalg.norm(Q-qa@(qa.T@Q))),passed=bool(ranks==[rank]*3 and max(np.linalg.norm(qa-Q@(Q.T@qa)),np.linalg.norm(Q-qa@(qa.T@Q)))<=1e-8))
def main():
    dest=ROOT/'rd/structure/main-shell/moment-enriched/records/moment-enriched-q0-space-v2-result.json';assert not dest.exists();records=[];stops=[]
    for n in [2,4]:
        try:
            start=time.monotonic();g=t.m.mesh(.4,n,.2);cells=[]
            for ids in g['cells']:
                q=g['ref'][ids].reshape(-1);op=t.operators(q);s=op['state'];Z,info=exact.local_basis(op,.4);nc=3*s['nscalar'];ni=3*s['ni']
                cells.append(dict(q=q,Z=Z,info=info,nc=nc,ni=ni,nb=nc-ni,ny=Z.shape[1],nscalar=s['nscalar'],faces=[dict(axis=x['axis'],sign=x['sign'],offset=x['offset'],n=x['n'],Mcoeff=x['Mcoeff']) for x in s['faces'][2:]]))
            Z,info,_,_=solver.assemble_constraints(cells,n);image=[];known=[];point_image=[];point_known=[]
            for cell in cells:
                ys=slice(cell['yoffset'],cell['yoffset']+cell['ny']);B=cell['Z'][:cell['nc']]@Z[ys];image.append(B);s=t.operators(cell['q'])['state'];v=t.f.reduce(s,lambda p:field(cell['q'],p,.4));known.append(v.reshape(s['nscalar'],6,3).transpose(0,2,1).reshape(cell['nc'],6))
                # Invert Ih on its full 36-dimensional q0 range, then inspect
                # the underlying polynomial at different points, not just traces.
                ids=[key for key in t.f.PID if key[2]==0];coeff=np.zeros((len(t.c.IDS),36))
                for j,key in enumerate(ids):coeff[t.c.LOOK[key],j]=1.
                I=np.zeros((s['nscalar'],36));I[:76]=coeff[t.f.VI]
                for face in s['faces']:I[face['offset']:face['offset']+face['n']]=face['Mcoeff'].T@t.c.restrict(coeff,face['axis'],face['sign'])
                modal=np.linalg.lstsq(I,B.reshape(s['nscalar'],-1),rcond=1e-12)[0];p,w=t.f.volume_points(15,9);L,_=t.c.scalar_basis(p,True);values=(L@coeff@modal).reshape(len(p),3,-1);point_image.append((values*np.sqrt(w/8)[:,None,None]).reshape(-1,Z.shape[1]))
                point_known.append((field(cell['q'],p,.4).reshape(len(p),6,3).transpose(0,2,1)*np.sqrt(w/8)[:,None,None]).reshape(-1,6))
            native=span(np.vstack(image),np.vstack(known),6);physical=span(np.vstack(point_image),np.vstack(point_known),6);seconds=time.monotonic()-start;peak=memory.memory_bytes()
            records.append(dict(n=n,q01_dimension=info['nullity'],native=native,polynomial=physical,seconds=seconds,peak_bytes=peak,passed=native['passed'] and physical['passed'] and seconds<=900 and peak<=2*1024**3));print(records[-1],flush=True)
        except Exception as exc:stops.append(dict(n=n,reason=repr(exc)));print('STOP',n,repr(exc),flush=True)
    paths=[Path(__file__),Path(solver.__file__),Path(exact.__file__),ROOT/'rd/structure/main-shell/moment-enriched/records/moment-enriched-thin-space-diagnostic-protocol.md',ROOT/'rd/structure/main-shell/moment-enriched/records/moment-enriched-thin-global-v2-result.json']
    dest.write_text(json.dumps(dict(records=records,stops=stops,six_field_hypothesis_verified=len(records)==2 and not stops and all(r['passed'] for r in records),geometry_conditional=True,finite_h_lower_bound=False,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}),indent=2),encoding='utf-8')
if __name__=='__main__':main()
