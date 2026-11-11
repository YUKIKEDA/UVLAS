"""Finite-h midsurface shear diagnostic; not a thin-limit locking theorem."""
from pathlib import Path
import hashlib,itertools,json,time
import numpy as np
from scipy.linalg import qr,lstsq
from numpy.polynomial.legendre import leggauss
import q9_edge_moment as method
import check_q9_curved_strain_screen as raw
import check_q9_distorted_plate as grid
import check_q9_independent_load as shape
import check_q9_plate_convergence as memory
ROOT=Path(__file__).resolve().parents[2];F=ROOT/'docs/verification';TH=[1e-9,1e-10,1e-11]
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def surface(name,x,y):
    dtype=np.result_type(x,y,float);r=np.array([x,y,0*x],dtype=dtype);rx=np.array([1+0*x,0*x,0*x],dtype=dtype);ry=np.array([0*x,1+0*x,0*x],dtype=dtype)
    if name=='SKEW':r[0]+=.4*y;ry[0]=.4
    elif name=='TRAP':r[0]*=1+.6*y;rx[0]=1+.6*y;ry[0]=.6*x
    elif name in ['SADDLE','NORMAL-C']:r[2]=.4*x*y;rx[2]=.4*y;ry[2]=.4*x
    normal=np.cross(rx,ry);normal=normal/np.sqrt(normal@normal)
    d=np.cross(rx,ry) if name=='SADDLE' else normal
    c=.2*d if name=='NORMAL-C' else np.zeros(3)
    return r,d,c,normal,rx,ry

def displacement(name,x,y):return np.sin(np.pi*x)*np.sin(np.pi*y)*surface(name,x,y)[3]

def make_mesh(name,n):
    vertices,xy,cells=grid.mesh(n,0. if name=='RECT' else .2);xy=xy-.5
    coords=[];normals=[];u=[];maximum=0.
    for x,y in xy:
        r,d,c,normal,rx,ry=surface(name,x,y)
        # The global parameter interval is [-.5,.5]; use shifted arguments for w.
        def value(a,b):return np.sin(np.pi*(a+.5))*np.sin(np.pi*(b+.5))*surface(name,a,b)[3]
        dx=value(x+1e-25j,y).imag/1e-25;dy=value(x,y+1e-25j).imag/1e-25
        A=np.array([rx,ry,d]);rhs=np.array([-d@dx,-d@dy,0.]);dd=np.linalg.solve(A,rhs)
        maximum=max(maximum,float(np.linalg.norm(A@dd-rhs)))
        coords.append([r,d,c]);normals.append(normal);u.append(np.array([value(x,y),dd,np.zeros(3)]))
    return xy,cells,np.array(coords),np.array(normals),np.array(u),maximum

def local(q,h,order,observation):
    state=method.prepare(q,h,0.,order);p,w=leggauss(observation)
    points=np.array([[x,y,0.] for x,y in itertools.product(p,p)]);weights=np.array([a*b for a,b in itertools.product(w,w)])
    J,E,B=method.evaluate(state,points);d=J[:,:,2]*2/h;a=J[:,:,:2];G=a.transpose(0,2,1)@a
    L=np.linalg.cholesky(G);area=np.sqrt(np.linalg.det(G));M=np.einsum('pia,kij,pj->pak',a,method.MB,d)*2
    white=np.linalg.solve(L,M)/np.linalg.norm(d,axis=1)[:,None,None]
    C=np.einsum('pak,pkn->pan',white,B[:,:,:81]);R=np.einsum('pak,pkn->pan',white,E[:,:,:81])
    factor=np.sqrt(weights*area)[:,None,None]
    return (C*factor).reshape(-1,81),(R*factor).reshape(-1,81),float(weights@area),state

def assemble(name,n,h,order=5,observation=5):
    xy,cells,coords,normals,u,analytic=make_mesh(name,n);side=2*n+1
    interior=np.array([j*side+i for j in range(1,side-1) for i in range(1,side-1)]);wm={v:i for i,v in enumerate(interior)}
    rows=2*observation**2;D=np.zeros((rows*n*n,3*len(xy)));W=np.zeros((len(D),len(interior)));known=[];original=[];area=0.;edge_ranks=set();edge_condition=0.
    for k,(_,_,ids,_) in enumerate(cells):
        C,R,A,state=local(coords[ids].reshape(81),h,order,observation);sl=slice(k*rows,(k+1)*rows);area+=A
        for j,node in enumerate(ids):
            D[sl,3*node:3*node+3]+=C[:,9*j+3:9*j+6]
            if node in wm:W[sl,wm[node]]+=C[:,9*j:9*j+3]@normals[node]/n
        known.append(C@u[ids].reshape(81));original.append(R@u[ids].reshape(81))
        for edge in state['edges'].values():edge_ranks.add(tuple(edge['ranks']));edge_condition=max(edge_condition,edge['condition'])
    D/=np.sqrt(area);W/=np.sqrt(area)
    return dict(D=D,W=W,known=np.concatenate(known)/np.sqrt(area),original=np.concatenate(original)/np.sqrt(area),xy=xy,interior=interior,u=u,normals=normals,analytic=analytic,edge_ranks=sorted(edge_ranks),edge_condition=edge_condition)

def main():
    target=F/'q9-edge-moment-shared-space-v2-result.json'
    if target.exists():raise FileExistsError(target)
    parent=F/'q9-edge-moment-field-audit.json';assert json.loads(parent.read_text())['all_criteria_pass']
    start=time.monotonic();out=ROOT/'build/analysis'/('q9-edge-moment-shared-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir();cases=[];stop=None;holdout=None
    try:
        for name,h,n in itertools.product(['RECT','SKEW','TRAP','SADDLE','NORMAL-C'],[.1,.01,.001],[1,2,4]):
            a=assemble(name,n,h);D,W=a['D'],a['W'];combined=np.column_stack((D,W));U,sd,_=np.linalg.svd(D,full_matrices=False);sc=np.linalg.svd(combined,compute_uv=False);scale=sc[0];ranks={}
            # Residual projection avoids constructing a square full left-null basis.
            for t in TH:
                rd=int(np.count_nonzero(sd>t*scale));rc=int(np.count_nonzero(sc>t*scale));Z=W-U[:,:rd]@(U[:,:rd].T@W);sz=np.linalg.svd(Z,compute_uv=False);rz=int(np.count_nonzero(sz>t*scale))
                ranks[str(t)]=dict(direction=rd,combined=rc,compatible_w=W.shape[1]-(rc-rd),projected_w=W.shape[1]-rz)
            qranks=[]
            for mat in [D,combined]:
                _,R,_=qr(mat,mode='economic',pivoting=True);qranks.append(int(np.count_nonzero(abs(R.diagonal())>1e-10*scale)))
            nodalw=np.sin(np.pi*(a['xy'][a['interior'],0]+.5))*np.sin(np.pi*(a['xy'][a['interior'],1]+.5))*n
            direction,_,_,_=lstsq(D,-W@nodalw,cond=1e-10,lapack_driver='gelsd');residual=D@direction+W@nodalw
            stable=len({tuple(r.values()) for r in ranks.values()})==1;central=ranks['1e-10']
            checks=dict(analytic=a['analytic']<=1e-11,rank_routes=all(r['compatible_w']==r['projected_w'] for r in ranks.values()),qr=qranks==[central['direction'],central['combined']],flat_w=name not in ['RECT','SKEW','TRAP'] or ranks['1e-09']['compatible_w']==W.shape[1],finite=bool(np.all(np.isfinite(direction))))
            path=out/f'{name}-h{h:g}-n{n}.npz';np.savez_compressed(path,D=D,W=W,known=a['known'],original=a['original'],xy=a['xy'],interior=a['interior'],direction_singular_values=sd,combined_singular_values=sc)
            cases.append(dict(geometry=name,h=h,n=n,w_dofs=W.shape[1],direction_dofs=D.shape[1],ranks=ranks,rank_stable=stable,qr_ranks=qranks,manufactured_norm=float(np.linalg.norm(a['known'])),raw_norm=float(np.linalg.norm(a['original'])),best_fit_norm=float(np.linalg.norm(residual)),best_direction_norm=float(np.linalg.norm(direction)),original_direction_norm=float(np.linalg.norm(a['u'][:,1,:])),edge_ranks=a['edge_ranks'],edge_condition=a['edge_condition'],checks={k:bool(v) for k,v in checks.items()},artifact=path.relative_to(ROOT).as_posix(),artifact_sha256=digest(path)))
            print(json.dumps({k:cases[-1][k] for k in ['geometry','h','n','rank_stable','manufactured_norm','best_fit_norm','checks']}),flush=True)
            if time.monotonic()-start>300 or memory.memory_bytes()>1024**3:raise RuntimeError('resource_limit')
        base=assemble('NORMAL-C',2,.001,5,5);fine=assemble('NORMAL-C',2,.001,7,5);obs=assemble('NORMAL-C',2,.001,5,9)
        opdiff=float(np.linalg.norm(np.column_stack((base['D']-fine['D'],base['W']-fine['W'])))/np.linalg.norm(np.column_stack((fine['D'],fine['W']))));normdiff=abs(np.linalg.norm(base['known'])/np.linalg.norm(obs['known'])-1)
        holdout=dict(construction_difference=opdiff,observation_norm_difference=float(normdiff),checks_pass=bool(opdiff<=1e-7 and normdiff<=1e-5))
        assert sum(p.stat().st_size for p in out.iterdir())<=32*1024**2
    except Exception as exc:stop=repr(exc)
    result=dict(source_sha256=digest(__file__),protocol_sha256=digest(F/'q9-edge-moment-shared-space-v2-experiment.md'),parent_sha256=digest(parent),dependencies={Path(m.__file__).name:digest(m.__file__) for m in [method,raw,grid,shape,memory]},cases=cases,holdout=holdout,stop_reason=stop,seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),diagnostic_checks_pass=bool(stop is None and len(cases)==45 and all(all(c['checks'].values()) for c in cases) and holdout['checks_pass']),product_shell_adopted=False)
    target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps({k:v for k,v in result.items() if k not in ['cases','dependencies']}))
    if not result['diagnostic_checks_pass']:raise SystemExit(1)
if __name__=='__main__':main()
