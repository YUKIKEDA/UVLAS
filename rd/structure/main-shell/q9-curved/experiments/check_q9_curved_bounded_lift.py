"""H1-bounded shared r/d/c lifts at fixed interpolated transverse displacement."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,itertools,json,time
import numpy as np
from scipy.linalg import svd,solve_triangular
from numpy.polynomial.legendre import leggauss
import q9_curved_constrained_gradient as method
import check_q9_edge_moment_shared_space_v2 as geometry
import check_q9_independent_load as shapes
import check_q9_plate_convergence as memory
from q9_edge_moment_fields import quadrature
ROOT=rd_paths.ROOT;F=rd_paths.RECORDS;LENGTH=1.;WAVES=[(1,1),(2,1),(1,2)]
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def fields(name,xy):
    r=np.zeros((len(xy),3,3));p=np.zeros((len(xy),6,3));error=0.
    for node,(x,y) in enumerate(xy):
        _,d,_,_,rx,ry=geometry.surface(name,x,y);H=np.array([rx,ry,d])
        for k,(a,b) in enumerate(WAVES):
            def value(x,y):return np.sin(a*np.pi*(x+.5))*np.sin(b*np.pi*(y+.5))*geometry.surface(name,x,y)[3]
            r[node,:,k]=value(x,y);dx=value(x+1e-25j,y).imag/1e-25;dy=value(x,y+1e-25j).imag/1e-25;rhs=np.array([-d@dx,-d@dy,0.]);dd=np.linalg.solve(H,rhs);p[node,:3,k]=dd;p[node,3:,k]=LENGTH*(.2 if name=='NORMAL-C' else 0.)*dd;error=max(error,float(np.linalg.norm(H@dd-rhs)))
    return r,p,error
def assemble(name,n,h,observation=3):
    xy,cells,coords,_,_,_=geometry.make_mesh(name,n);r,reference,analytic=fields(name,xy);points,w=quadrature(observation);nr=2*len(points);A=np.zeros((nr*n*n,6*len(xy)));force=np.zeros((len(A),3));rawknown=np.zeros_like(force);M=np.zeros((len(xy),len(xy)));volume=0.;area=0.;edge_condition=0.
    p,gw=leggauss(5)
    for cell,(_,_,ids,_) in enumerate(cells):
        q=coords[ids].reshape(81);state=method.prepare(q,h,0.,9);J,E,B=method.evaluate(state,points);a=J[:,:,:2];dz=J[:,:,2]*2/h;G=a.transpose(0,2,1)@a;LG=np.linalg.cholesky(G);row=2*np.einsum('pia,kij,pj->pak',a,method.MB,dz);white=np.linalg.solve(LG,row)/np.linalg.norm(dz,axis=1)[:,None,None];weight=w*np.linalg.det(J);volume+=weight.sum();C=np.einsum('pak,pkn->pan',white,B[:,:,:81]);RC=np.einsum('pak,pkn->pan',white,E[:,:,:81]);C=(C*np.sqrt(weight)[:,None,None]).reshape(nr,81);RC=(RC*np.sqrt(weight)[:,None,None]).reshape(nr,81);sl=slice(nr*cell,nr*(cell+1))
        ur=np.zeros((9,9,3));ur[:,:3]=r[ids];force[sl]=C@ur.reshape(81,3);ur[:,3:6]=reference[ids,:3];ur[:,6:9]=reference[ids,3:]/LENGTH;rawknown[sl]=RC@ur.reshape(81,3)
        for local,node in enumerate(ids):A[sl,6*node:6*node+3]+=C[:,9*local+3:9*local+6];A[sl,6*node+3:6*node+6]+=C[:,9*local+6:9*local+9]/LENGTH
        localM=np.zeros((9,9))
        for i,j in itertools.product(range(5),repeat=2):
            N,dx,dy=shapes.shape(p[i],p[j]);g=np.column_stack((dx,dy));surfaceJ=coords[ids,0,:].T@g;metric=surfaceJ.T@surfaceJ;aw=gw[i]*gw[j]*np.sqrt(np.linalg.det(metric));area+=aw;localM+=aw*(np.outer(N,N)+LENGTH**2*g@np.linalg.inv(metric)@g.T)
        M[np.ix_(ids,ids)]+=localM
        edge_condition=max(edge_condition,max(e['condition'] for e in state['edges'].values()))
    A/=np.sqrt(volume);force/=np.sqrt(volume);rawknown/=np.sqrt(volume);M/=area
    return dict(A=A,F=force,M=M,xy=xy,reference=reference,analytic=analytic,edge_condition=edge_condition,raw_known=rawknown)
def whiten(A,M):
    chol=np.linalg.cholesky(M);nr=A.shape[0];N=len(M);AT=A.reshape(nr,N,6).transpose(1,0,2).reshape(N,-1);temp=solve_triangular(chol,AT,lower=True,check_finite=False);return temp.reshape(N,nr,6).transpose(1,0,2).reshape(nr,6*N),chol
def physical(v,chol):return solve_triangular(chol.T,v.reshape(len(chol),6),lower=False,check_finite=False).reshape(-1)
def pnorm(p,M):return float(np.sqrt(np.einsum('nc,nm,mc->',p,M,p)))
def bounded(U,s,Vh,f,bound,Aw,chol):
    rhs=U.T@f;mask=s>np.finfo(float).eps*max(Aw.shape)*s[0];coef=np.zeros_like(s);coef[mask]=-rhs[mask]/s[mask];lam=0.
    if np.linalg.norm(coef)>bound:
        low=0.;high=s[0]**2
        def coefficients(x):return -s*rhs/(s*s+x)
        while np.linalg.norm(coefficients(high))>bound:high*=4
        for _ in range(90):
            mid=(low+high)/2
            if np.linalg.norm(coefficients(mid))>bound:low=mid
            else:high=mid
        lam=high;coef=coefficients(lam)
    v=Vh.T@coef;residual=Aw@v+f;kkt=float(np.linalg.norm(Aw.T@residual+lam*v)/max(s[0]*np.linalg.norm(f)+lam*np.linalg.norm(v),1e-30));error=float(max(0.,np.linalg.norm(v)/bound-1))
    return dict(p=physical(v,chol),residual=float(np.linalg.norm(residual)),h1=float(np.linalg.norm(v)),direction_h1=float(np.linalg.norm(v.reshape(len(chol),6)[:,:3])),curvature_h1=float(np.linalg.norm(v.reshape(len(chol),6)[:,3:])),multiplier=float(lam),kkt=kkt,constraint_error=error)
def main():
    target=F/'q9-curved-bounded-lift-result.json'
    if target.exists():raise FileExistsError(target)
    parent=F/'q9-curved-constrained-gradient-audit.json';assert json.loads(parent.read_text())['all_criteria_pass'];start=time.monotonic();out=ROOT/'build/analysis'/('q9-curved-bounded-lift-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir();cases=[];holdout=None;stop=None;selected=None
    configs=[(name,n,h) for name,n,h in itertools.product(['SADDLE','NORMAL-C'],[2,4],[.1,.01,.001])]+[(name,8,.001) for name in ['SADDLE','NORMAL-C']]+[(name,n,h) for name,n,h in itertools.product(['RECT','TRAP'],[2,4],[.01,.001])]
    try:
        for name,n,h in configs:
            a=assemble(name,n,h);Aw,chol=whiten(a['A'],a['M']);U,s,Vh=svd(Aw,full_matrices=False,check_finite=False,lapack_driver='gesdd');rows=[];saved_bounded=[];saved_free=[]
            for k,wave in enumerate(WAVES):
                f=a['F'][:,k];reference=pnorm(a['reference'][:,:,k],a['M']);forcing=float(np.linalg.norm(f));free=[];fits=[];rhs=U.T@f
                for threshold in [1e-9,1e-10,1e-11]:
                    keep=s>threshold*s[0];coef=np.zeros_like(s);coef[keep]=-rhs[keep]/s[keep];v=Vh.T@coef;res=float(np.linalg.norm(Aw@v+f));free.append(dict(threshold=threshold,rank=int(keep.sum()),residual=res,relative=res/forcing,h1=float(np.linalg.norm(v)),direction_h1=float(np.linalg.norm(v.reshape(len(chol),6)[:,:3])),curvature_h1=float(np.linalg.norm(v.reshape(len(chol),6)[:,3:]))));saved_free.append(physical(v,chol))
                for multiplier in [2.,4.]:
                    fit=bounded(U,s,Vh,f,multiplier*reference,Aw,chol);saved_bounded.append(fit.pop('p'));fit.update(bound_factor=multiplier,bound=multiplier*reference,relative=fit['residual']/forcing,thin_scaled=fit['residual']/forcing*LENGTH/h);fits.append(fit)
                pref=a['reference'][:,:,k].reshape(-1);known=float(np.linalg.norm(a['A']@pref+f));rows.append(dict(wave=list(wave),forcing=forcing,reference_h1=reference,reference_residual=known,raw_reference_residual=float(np.linalg.norm(a['raw_known'][:,k])),unbounded=free,bounded=fits))
            path=out/f'{name}-n{n}-h{h:g}.npz';np.savez_compressed(path,A=a['A'],F=a['F'],M=a['M'],xy=a['xy'],reference=a['reference'],bounded=np.array(saved_bounded),unbounded=np.array(saved_free),singular_values=s)
            checks=dict(analytic=a['analytic']<=1e-11,finite=bool(np.all(np.isfinite(saved_bounded))),constraint=all(fit['constraint_error']<=1e-10 for row in rows for fit in row['bounded']),kkt=all(fit['kkt']<=1e-9 for row in rows for fit in row['bounded']))
            cases.append(dict(geometry=name,n=n,h=h,unknowns=Aw.shape[1],observations=Aw.shape[0],edge_condition=a['edge_condition'],analytic=a['analytic'],rows=rows,checks={k:bool(v) for k,v in checks.items()},artifact=path.relative_to(ROOT).as_posix(),artifact_sha256=digest(path)))
            if name=='NORMAL-C' and n==2 and h==.001:selected=(np.array(saved_bounded),a['M'].copy(),[fit['residual'] for row in rows for fit in row['bounded']])
            print(json.dumps(dict(geometry=name,n=n,h=h,checks=cases[-1]['checks'],bounded_relative=[row['bounded'][0]['relative'] for row in rows],free_h1=[row['unbounded'][1]['h1'] for row in rows])),flush=True)
            del Aw,U,Vh
            if time.monotonic()-start>600 or memory.memory_bytes()>1024**3:raise RuntimeError('resource_limit')
        other=assemble('NORMAL-C',2,.001,observation=5);solutions,metric,previous=selected;differences=[]
        for i,p in enumerate(solutions):
            value=float(np.linalg.norm(other['A']@p+other['F'][:,i//2]));differences.append(abs(value-previous[i])/max(value,previous[i],1e-30))
        holdout=dict(observation=5,relative_differences=differences,maximum=float(max(differences)),checks_pass=bool(max(differences)<=1e-4),metric_difference=float(np.linalg.norm(metric-other['M'])/np.linalg.norm(metric)))
        assert sum(p.stat().st_size for p in out.iterdir())<=64*1024**2
    except Exception as exc:stop=repr(exc)
    result=dict(source_sha256=digest(__file__),parent_sha256=digest(parent),contract_sha256=digest(ROOT/'rd/structure/main-shell/q9-curved/q9-curved-bounded-lift-contract.md'),protocol_sha256=digest(F/'q9-curved-bounded-lift-experiment.md'),dependencies={Path(m.__file__).name:digest(m.__file__) for m in [method,method.base,geometry,geometry.grid,shapes,memory]},cases=cases,holdout=holdout,stop_reason=stop,seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),diagnostic_checks_pass=bool(stop is None and len(cases)==22 and all(all(c['checks'].values()) for c in cases) and holdout['checks_pass']),product_shell_adopted=False)
    target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps({k:v for k,v in result.items() if k not in ['cases','dependencies']}))
    if not result['diagnostic_checks_pass']:raise SystemExit(1)
if __name__=='__main__':main()
