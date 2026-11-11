"""Reference-only coefficient audit, with protocol fixed before execution."""
from pathlib import Path
from decimal import Decimal, localcontext
import datetime
import gzip
import hashlib
import json
import math
import time
import numpy as np
from scipy.integrate import quad_vec
from scipy.linalg import solve
import check_avl_collinear_kernel as independent
import check_q9_plate_convergence as memory

ROOT=Path(__file__).resolve().parents[2];FOLDER=ROOT/'docs/verification'
PI=independent.PI
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dec(v):return Decimal.from_float(float(v))
def dot3(a,b):return (a[...,0]*b[...,0]+a[...,1]*b[...,1])+a[...,2]*b[...,2]
def norm(a):return np.sqrt(dot3(a,a))

def native_coefficients(p,a,b,normal):
    """Mach 0, zero core, all bound legs included: current VORVELC path."""
    aa=a-p;bb=b-p;asq=dot3(aa,aa);bsq=dot3(bb,bb);ra=np.sqrt(asq);rb=np.sqrt(bsq)
    edge=b-a;length=norm(edge);unit=edge/length[:,None]
    q=-dot3(aa,unit);outside=(q<0)|(q>length)
    perp=-aa-q[:,None]*unit
    bound=np.zeros_like(aa)
    k=outside
    si=length[k]*(2*q[k]-length[k])/(ra[k]*rb[k]*(q[k]*rb[k]+(q[k]-length[k])*ra[k]))
    bound[k]=np.cross(unit[k],perp[k])*si[:,None]
    k=~outside;cross=np.cross(aa[k],bb[k]);csq=dot3(cross,cross);ab=dot3(aa[k],bb[k])
    if np.any(csq==0):raise ValueError('contact_or_collinear_interior')
    factor=((bsq[k]-ab)/np.sqrt(np.sqrt(bsq[k]**2))+(asq[k]-ab)/np.sqrt(np.sqrt(asq[k]**2)))/np.sqrt(csq**2)
    bound[k]=cross*factor[:,None]
    for r,radius,sign in [(aa,ra,-1.),(bb,rb,1.)]:
        rsq=r[:,2]**2+r[:,1]**2
        if np.any(rsq==0) or np.any(radius==0):raise ValueError('trailing_line_contact')
        t=sign*(1-r[:,0]/radius)/np.sqrt(rsq**2)
        bound[:,1]+=r[:,2]*t;bound[:,2]-=r[:,1]*t
    coefficient=dot3(bound*PI,normal)
    finite_distance=norm(p-(a+np.clip(q,0,length)[:,None]*unit))
    distances=[finite_distance]
    for endpoint in [a,b]:
        r=p-endpoint;distances.append(np.sqrt(r[:,1]**2+r[:,2]**2+np.minimum(r[:,0],0)**2))
    return coefficient,np.minimum.reduce(distances)/length

def coefficient(p,a,b,normal,precision):
    vector=independent.decimal_pair(p,a,b,include=True,precision=precision)
    with localcontext() as context:
        context.prec=precision
        return +sum((v*dec(n) for v,n in zip(vector,normal)),Decimal(0))

def line_integral(p,a,b,normal):
    edge=b-a
    def integrand(t):
        r=p-(a+t*edge)
        return np.cross(edge,r)/np.linalg.norm(r)**3
    finite,error,info=quad_vec(integrand,0.,1.,epsabs=1e-12,epsrel=1e-12,full_output=True)
    total=finite.copy();errors=[float(error)];ok=bool(info.success);evaluations=int(info.neval)
    for endpoint,sign in [(a,-1.),(b,1.)]:
        r=p-endpoint;rho=float(np.linalg.norm(r[1:]))
        if rho==0:raise ValueError('trailing_line_contact')
        lower=math.atan2(-r[0],rho)
        direction=sign*np.cross(np.array([1.,0.,0.]),r)/rho**2
        value,error,info=quad_vec(lambda theta:direction*math.cos(theta),lower,math.pi/2,
            epsabs=1e-12,epsrel=1e-12,full_output=True)
        total+=value;errors.append(float(error));ok=ok and bool(info.success);evaluations+=int(info.neval)
    return dict(coefficient=float(np.dot(normal,PI*total)),estimated_velocity_error=PI*math.fsum(errors),
        integrator_success=ok,evaluations=evaluations)

def small_solve(a80,original,rows,vectors,adjoints):
    b=-rows[:,6];base=solve(original,b);candidate=solve(a80,b)
    w=adjoints[:,:,2]
    project=lambda g:np.array([math.fsum(float(x)*float(y) for x,y in zip(weight,g)) for weight in w])
    qbase=project(base);qnew=project(candidate);qstored=project(vectors[:,1]);qir=project(vectors[:,4])
    metrics=[]
    for q in [qbase,qstored,qir]:
        d=qnew-q
        metrics.append(dict(lift=float(abs(d[0])/abs(qstored[0])),moment_y=float(abs(d[1])/max(abs(qstored[1]),.1)),
            lift_bands=float(np.linalg.norm(d[2:])/np.linalg.norm(qstored[2:]))))
    predicted=-adjoints[:,:,3]@((a80-original)@base)
    return dict(operator='Decimal80 coefficients projected to float64; SciPy float64 solve',
        q_original_scipy=qbase.tolist(),q_reference_scipy=qnew.tolist(),q_previous_original=qstored.tolist(),
        q_previous_ir2=qir.tolist(),relative_changes=metrics,coefficient_first_order_prediction=predicted.tolist(),
        actual_minus_prediction=(qnew-qbase-predicted).tolist(),
        candidate_residual_max=float(np.max(abs(b-a80@candidate))),
        original_scipy_residual_max=float(np.max(abs(b-original@base))))

def analyze(run,root,begin):
    folder=ROOT/run['directory']/'out';rows=np.loadtxt(folder/'uv-rows.txt',skiprows=2)
    vectors=np.loadtxt(folder/'us-vectors.txt',skiprows=1);n=len(rows)
    adj=np.loadtxt(folder/'us-adjoints.txt',skiprows=1).reshape((22,n,5))
    sel=np.array([0,1,n//4-1,n//2-1,n//2,3*n//4-1,n-2,n-1])
    if n<=192:
        original=np.fromfile(folder/'us-matrix.bin',dtype='<f8').reshape((n,n),order='F')
        ii,jj=np.indices((n,n));ii=ii.ravel();jj=jj.ravel();stored=original[ii,jj]
        reference_matrix=np.empty_like(original)
    else:
        cross=np.loadtxt(folder/'us-matrix-cross.txt',skiprows=1)
        ii=cross[:,0].astype(int)-1;jj=cross[:,1].astype(int)-1;stored=cross[:,2]
        reference_matrix=None
    assert len(set(zip(ii.tolist(),jj.tolist())))==len(ii)
    assert np.all(np.isfinite(rows)) and np.max(abs(rows[:,4]))==0
    points=rows[ii,1:4];a=rows[jj,7:10];b=rows[jj,10:13];normals=rows[ii,4:7]
    native,distance=native_coefficients(points,a,b,normals)
    out=root/f'{run["family"]}-{run["case"]}-M{run["m"]}-N{run["n"]}';out.mkdir()
    coefficient_file=out/'coefficients.tsv.gz';differences=np.empty(len(ii));projected=np.empty(len(ii))
    g=[dec(v) for v in vectors[:,4]];row_action=[Decimal(0) for _ in range(n)];row_absolute=row_action.copy()
    with localcontext() as context,gzip.open(coefficient_file,'wt',encoding='utf-8',newline='\n') as stream:
        context.prec=80
        stream.write('i\tj\tstored_float64\treference_decimal80\tdelta_reference_minus_stored\n')
        for k,(i,j) in enumerate(zip(ii,jj)):
            value=coefficient(points[k],a[k],b[k],normals[k],80);delta=value-dec(stored[k])
            differences[k]=float(delta);projected[k]=float(value)
            term=delta*g[j];row_action[i]+=term;row_absolute[i]+=abs(term)
            stream.write(f'{i+1}\t{j+1}\t{stored[k]:.17e}\t{value}\t{delta}\n')
            if k%4096==0 and time.monotonic()-begin>300:raise TimeoutError('coefficient_audit_time_limit')
    target=1e-10+1e-10*abs(projected);ratio=abs(differences)/target
    order=sorted(range(len(ii)),key=lambda k:(-ratio[k],int(ii[k]),int(jj[k])))
    nearest=sorted(range(len(ii)),key=lambda k:(distance[k],int(ii[k]),int(jj[k])))[:16]
    selected=sorted(set(order[:16]+nearest));selfaudit=[]
    for k in selected:
        low=coefficient(points[k],a[k],b[k],normals[k],80)
        high=coefficient(points[k],a[k],b[k],normals[k],160)
        with localcontext() as ctx:
            ctx.prec=160;difference=abs(high-low)
            eligible=difference<=Decimal('1e-20') or difference<=Decimal('1e-20')*abs(high)
        selfaudit.append(dict(i=int(ii[k]+1),j=int(jj[k]+1),low=str(low),high=str(high),difference=str(difference),eligible=eligible))
    diagonal=int(np.flatnonzero((ii==0)&(jj==0))[0]);other=next(k for k in order if k!=diagonal)
    quadrature=[]
    for k in [diagonal,other]:
        value=line_integral(points[k],a[k],b[k],normals[k])
        value.update(i=int(ii[k]+1),j=int(jj[k]+1),reference=float(projected[k]),difference=abs(value['coefficient']-projected[k]))
        value['eligible']=bool(value['integrator_success'] and value['difference']<=target[k])
        quadrature.append(value)
    selected_rows=[]
    for i in sel:
        selected_rows.append(dict(i=int(i+1),coefficient_defect_action=str(row_action[i]),
            absolute_term_sum=str(row_absolute[i]),previous_stored_residual=float(vectors[i,8])))
    with localcontext() as ctx:
        ctx.prec=80
        partial=[-sum((dec(z)*v for z,v in zip(adj[k,:,3],row_action)),Decimal(0)) for k in range(22)]
    worst=[]
    for k in order[:16]:
        worst.append(dict(i=int(ii[k]+1),j=int(jj[k]+1),stored=float(stored[k]),reference=float(projected[k]),
            difference=float(differences[k]),scaled_error=float(ratio[k]),relative_geometric_distance=float(distance[k])))
    small=None
    if reference_matrix is not None:
        reference_matrix[ii,jj]=projected
        small=small_solve(reference_matrix,original,rows,vectors,adj)
        np.savez_compressed(out/'projected-small-matrix.npz',reference=reference_matrix)
    np.savez_compressed(out/'coefficient-diagnostics.npz',i=ii+1,j=jj+1,stored=stored,
        reference_projected=projected,decimal_difference=differences,reproduced=native,relative_geometric_distance=distance)
    records={p.name:dict(sha256=digest(p),bytes=p.stat().st_size) for p in out.iterdir()}
    rec=dict(family=run['family'],case=run['case'],m=run['m'],n=run['n'],unknowns=n,
        directory=out.relative_to(ROOT).as_posix(),coefficient_count=len(ii),full_matrix=n<=192,
        maximum_coefficient_difference=float(max(abs(differences))),maximum_scaled_error=float(max(ratio)),
        maximum_native_reproduction_difference=float(max(abs(native-stored))),
        native_reproduction_eligible=bool(np.all(abs(native-stored)<=1e-10+1e-10*abs(stored))),
        coefficient_eligible=bool(np.all(abs(differences)<=target)),minimum_relative_distance=float(min(distance)),
        diagonal_includes_bound=True,decimal_self_audits=selfaudit,independent_line_integrals=quadrature,
        selected_rows=selected_rows,partial_first_order_output_defect=[str(v) for v in partial],
        partial_is_full=n<=192,partial_is_error_bound=False,worst_pairs=worst,small_solve=small,files=records)
    summary=out/'case-summary.json'
    summary.write_text(json.dumps(rec,indent=2,allow_nan=False)+'\n')
    records[summary.name]=dict(sha256=digest(summary),bytes=summary.stat().st_size)
    print(json.dumps(dict(family=run['family'],case=run['case'],m=run['m'],n=run['n'],pairs=len(ii),
        max_error=rec['maximum_coefficient_difference'],eligible=rec['coefficient_eligible'],seconds=time.monotonic()-begin)),flush=True)
    return rec

def main():
    begin=time.monotonic();path=FOLDER/'avl-solve-sensitivity-runs.json';parent=json.loads(path.read_text())
    output=FOLDER/'avl-aic-generation-result.json'
    if output.exists():raise FileExistsError(output)
    root=ROOT/'build/analysis'/('avl-aic-generation-'+datetime.datetime.now().strftime('%Y%m%d-%H%M%S'));root.mkdir()
    records=[];stop=None
    try:
        for run in parent['runs']:
            records.append(analyze(run,root,begin))
            (root/'progress.json').write_text(json.dumps(dict(completed=len(records)))+'\n')
            if sum(p.stat().st_size for p in root.rglob('*') if p.is_file())>128*1024**2:raise RuntimeError('disk_budget')
            if memory.memory_bytes()>512*1024**2:raise RuntimeError('memory_budget')
        checks=dict(eight_systems=len(records)==8,
            actual_coefficients_reproduced=all(r['native_reproduction_eligible'] for r in records),
            decimal_coefficient_accuracy=all(r['coefficient_eligible'] for r in records),
            decimal_self_audits=all(v['eligible'] for r in records for v in r['decimal_self_audits']),
            independent_line_integrals=all(v['eligible'] for r in records for v in r['independent_line_integrals']),
            full_small_matrix_output_precision=all(max(max(v.values()) for v in r['small_solve']['relative_changes'])<=1e-10 for r in records if r['full_matrix']),
            resource_budget=time.monotonic()-begin<300 and memory.memory_bytes()<512*1024**2 and sum(p.stat().st_size for p in root.rglob('*') if p.is_file())<128*1024**2)
    except Exception as exc:
        stop=repr(exc);checks=dict(completed=False);raise
    finally:
        result=dict(protocol='docs/verification/avl-aic-generation-experiment.md',source_sha256=digest(Path(__file__)),
            dependency_sha256={Path(m.__file__).name:digest(Path(m.__file__)) for m in [independent,memory]},
            parent_runs_sha256=digest(path),protocol_sha256=digest(FOLDER/'avl-aic-generation-experiment.md'),
            kernel_sha256=digest(ROOT/'build/reference-avl/source-core0-grid/src/aic.f'),
            setup_sha256=digest(ROOT/'build/reference-avl/source-core0-grid/src/asetup.f'),
            records=records,checks=checks,all_criteria_pass=all(checks.values()),stop_reason=stop,
            reference_convergence_accepted=False,large_matrix_global_error_bounded=False,
            seconds=time.monotonic()-begin,peak_process_working_bytes=memory.memory_bytes(),
            total_bytes=sum(p.stat().st_size for p in root.rglob('*') if p.is_file()))
        output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(checks=checks,seconds=result['seconds'],peak=result['peak_process_working_bytes'],bytes=result['total_bytes']),indent=2))

if __name__=='__main__':main()
