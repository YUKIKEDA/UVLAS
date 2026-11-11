"""Predeclared native-pair reproduction and Decimal collinear-kernel audit."""
from pathlib import Path
from decimal import Decimal, localcontext
import hashlib
import json
import math
import time
import numpy as np
import check_q9_plate_convergence as memory

ROOT=Path(__file__).resolve().parents[2]
START=time.monotonic()
PI=.079577472
RECORDS=[]
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def norm(x):return np.linalg.norm(x,axis=-1)

def native_parts(points,a,b,exclude):
    aa=a[None,:,:]-points[:,None,:];bb=b[None,:,:]-points[:,None,:]
    asq=np.sum(aa*aa,axis=-1);bsq=np.sum(bb*bb,axis=-1)
    ar=np.sqrt(asq);br=np.sqrt(bsq);cross=np.cross(aa,bb)
    csq=np.sum(cross*cross,axis=-1);dot=np.sum(aa*bb,axis=-1)
    numerator=(bsq-dot)/np.sqrt(np.sqrt(bsq**2))+(asq-dot)/np.sqrt(np.sqrt(asq**2))
    factor=np.divide(numerator,np.sqrt(csq**2),out=np.zeros_like(csq),where=(csq!=0)&~exclude)
    bound=cross*factor[:,:,None]
    ta=-(1-aa[:,:,0]/ar)/np.sqrt(np.sum(aa[:,:,1:]**2,axis=-1)**2)
    tb=(1-bb[:,:,0]/br)/np.sqrt(np.sum(bb[:,:,1:]**2,axis=-1)**2)
    trail=np.zeros_like(bound)
    trail[:,:,1]=aa[:,:,2]*ta+bb[:,:,2]*tb
    trail[:,:,2]=-aa[:,:,1]*ta-bb[:,:,1]*tb
    return PI*bound,PI*trail

def stable_bound(points,a,b,exclude):
    edge=b-a;length=norm(edge);unit=edge/length[:,None]
    r=points[:,None,:]-a[None,:,:];q=np.sum(r*unit[None,:,:],axis=-1)
    perp=r-q[:,:,None]*unit[None,:,:];d2=np.sum(perp*perp,axis=-1)
    ra=norm(r);rb=norm(points[:,None,:]-b[None,:,:]);v=q-length[None,:]
    outside=(q<0)|(v>0)
    if np.any((~outside)&(d2==0)&~exclude):raise ValueError('finite_segment_contact')
    integ=np.zeros_like(q)
    denominator=ra*rb*(q*rb+v*ra)
    np.divide(length[None,:]*(2*q-length[None,:]),denominator,out=integ,where=outside&~exclude)
    interior=~outside&~exclude
    np.divide(q/ra-v/rb,d2,out=integ,where=interior)
    return PI*np.cross(unit[None,:,:],perp)*integ[:,:,None],dict(q=q,length=length,distance=np.sqrt(d2),outside=outside)

def stable_trail(points,a,b):
    result=np.zeros((len(points),len(a),3))
    for endpoint,sign in [(a,-1.),(b,1.)]:
        r=endpoint[None,:,:]-points[:,None,:];radius=norm(r);d2=np.sum(r[:,:,1:]**2,axis=-1)
        upstream=r[:,:,0]>=0
        t=np.empty_like(radius)
        np.divide(sign,radius*(radius+r[:,:,0]),out=t,where=upstream)
        np.divide(sign*(1-r[:,:,0]/radius),d2,out=t,where=~upstream)
        result[:,:,1]+=r[:,:,2]*t
        result[:,:,2]-=r[:,:,1]*t
    return PI*result

def decimal_pair(p,a,b,include=True,precision=80,bound_only=False):
    with localcontext() as ctx:
        ctx.prec=precision
        cast=lambda v:Decimal.from_float(float(v))
        aa=[cast(x)-cast(y) for x,y in zip(a,p)]
        bb=[cast(x)-cast(y) for x,y in zip(b,p)]
        dot=lambda x,y:sum((u*v for u,v in zip(x,y)),Decimal(0))
        asq=dot(aa,aa);bsq=dot(bb,bb);ar=asq.sqrt();br=bsq.sqrt()
        cross=[aa[1]*bb[2]-aa[2]*bb[1],aa[2]*bb[0]-aa[0]*bb[2],aa[0]*bb[1]-aa[1]*bb[0]]
        csq=dot(cross,cross);ab=dot(aa,bb);out=[Decimal(0)]*3
        if include and csq:
            factor=((bsq-ab)/br+(asq-ab)/ar)/csq
            out=[c*factor for c in cross]
        elif include and not csq:
            edge=[y-x for x,y in zip(aa,bb)];q=-dot(aa,edge)/dot(edge,edge)
            if 0<=q<=1:raise ValueError('finite_segment_contact')
        if not bound_only:
            for r,radius,sign in [(aa,ar,-1),(bb,br,1)]:
                d2=r[1]**2+r[2]**2
                t=Decimal(sign)*(1-r[0]/radius)/d2
                out[1]+=r[2]*t;out[2]-=r[1]*t
        return [+(v*cast(PI)) for v in out]

def weighted_sum(kernel,gamma):
    return np.array([[math.fsum(float(v)*float(g) for v,g in zip(row[:,k],gamma)) for k in range(3)] for row in kernel])

def force_coefficients(forces,bound,origin,s,c,b):
    changed=forces.copy();mask=changed[:,2]==0
    changed[mask,6:9]=bound[changed[mask,0].astype(int)-1]
    f=np.array([math.fsum(changed[:,6+j]) for j in range(3)])/s
    moments=np.cross(changed[:,3:6]-origin,changed[:,6:9])
    m=np.array([math.fsum(moments[:,j]) for j in range(3)])/(s*np.array([b,c,b]))
    return dict(force=f.tolist(),moment=m.tolist())

def main():
    parent_path=ROOT/'docs/verification/avl-pair-probe-runs.json';parent=json.loads(parent_path.read_text())
    for run in parent['runs']:
        folder=ROOT/run['directory'];m=run['m'];n=run['n'];size=m*n
        points_data=np.loadtxt(folder/'uv-force-points.txt',skiprows=2)
        rows=np.loadtxt(folder/'uv-rows.txt',skiprows=2)
        forces=np.loadtxt(folder/'uv-forces.txt',skiprows=2)
        native_dump=np.loadtxt(folder/'uv-force-matrix.txt',skiprows=1)[:,2:].reshape(size,size,3)
        points=points_data[:,1:4];gamma=points_data[:,4];known=points_data[:,5:8];vv=points_data[:,8:11]
        a=rows[:,7:10];b=rows[:,10:13];exclude=np.eye(size,dtype=bool)
        old_bound,old_trail=native_parts(points,a,b,exclude)
        new_bound,geometry=stable_bound(points,a,b,exclude);new_trail=stable_trail(points,a,b)
        old=old_bound+old_trail;new=new_bound+new_trail
        old_velocity=weighted_sum(old,gamma);stable_velocity=weighted_sum(new,gamma)
        reference=np.empty_like(old)
        for i in range(size):
            for j in range(size):
                reference[i,j]=[float(v) for v in decimal_pair(points[i],a[j],b[j],i!=j)]
            if time.monotonic()-START>180:raise TimeoutError('decimal_trial_time_limit')
        actual=norm(new-reference);target=1e-10+1e-10*norm(reference)
        discrepancy=norm(old-reference)
        worst_indices=np.argsort(discrepancy.reshape(-1))[-16:][::-1]
        worst=[];decimal_checks=[]
        for index in worst_indices:
            i,j=np.unravel_index(index,(size,size))
            lo=decimal_pair(points[i],a[j],b[j],i!=j,80)
            hi=decimal_pair(points[i],a[j],b[j],i!=j,160)
            with localcontext() as context:
                context.prec=160
                err=sum(((x-y)**2 for x,y in zip(lo,hi)),Decimal(0)).sqrt()
                scale=sum((x*x for x in hi),Decimal(0)).sqrt()
                eligible=err<=Decimal('1e-20') or err<=Decimal('1e-20')*scale
            decimal_checks.append(dict(i=int(i+1),j=int(j+1),low=[str(x) for x in lo],high=[str(x) for x in hi],difference=str(err),eligible=eligible))
            worst.append(dict(i=int(i+1),j=int(j+1),source_a=a[j].tolist(),source_b=b[j].tolist(),point=points[i].tolist(),
                projection_fraction=float(geometry['q'][i,j]/geometry['length'][j]),
                relative_distance=float(geometry['distance'][i,j]/geometry['length'][j]),outside=bool(geometry['outside'][i,j]),
                original=old[i,j].tolist(),stable=new[i,j].tolist(),reference=reference[i,j].tolist(),
                original_error=float(discrepancy[i,j]),stable_error=float(actual[i,j])))
        old_force=2*gamma[:,None]*np.cross(known+old_velocity,b-a)
        native_bound_force=forces[forces[:,2]==0,6:9]
        stable_force=2*gamma[:,None]*np.cross(known+stable_velocity,b-a)
        bound_only_velocity=weighted_sum(new_bound+old_trail,gamma)
        bound_only_force=2*gamma[:,None]*np.cross(known+bound_only_velocity,b-a)
        origin=np.array([.375,0.,0.]);s=8.4;c=1.1142857142857143;span=8.
        old_coeff=force_coefficients(forces,old_force,origin,s,c,span)
        new_coeff=force_coefficients(forces,stable_force,origin,s,c,span)
        bound_coeff=force_coefficients(forces,bound_only_force,origin,s,c,span)
        mirror=np.array([(n-1-j)*m+k for j in range(n) for k in range(m)])
        reflection=np.array([1.,-1.,1.])
        geom_error=max(np.max(norm(points-points[mirror]*reflection)),np.max(norm(a-b[mirror]*reflection)),np.max(norm(b-a[mirror]*reflection)))
        gamma_error=float(np.max(abs(gamma-gamma[mirror])/np.maximum(abs(gamma),1.)))
        record=dict(case=run['case'],m=m,n=n,pair_count=size**2,
            native_matrix_reproduction=float(np.max(norm(old-native_dump))),
            native_velocity_reproduction=float(np.max(norm(old_velocity-vv))),
            native_force_reproduction=float(np.max(abs(old_force-native_bound_force))/s),
            stable_accuracy_eligible=bool(np.all(actual<=target)),maximum_stable_error=float(np.max(actual)),
            maximum_native_error=float(np.max(discrepancy)),decimal_comparisons=decimal_checks,worst_pairs=worst,
            original_coefficients=old_coeff,stable_coefficients=new_coeff,bound_only_coefficients=bound_coeff,
            geometry_mirror_error=float(geom_error),gamma_mirror_error=gamma_error,
            bound_replacement_explains_lateral=abs(bound_coeff['force'][1])<=1e-9 and abs(bound_coeff['moment'][2])<=1e-9)
        RECORDS.append(record)
        print(json.dumps(dict(case=run['case'],m=m,n=n,native_reproduction=record['native_velocity_reproduction'],stable_error=record['maximum_stable_error'],stable_CY=new_coeff['force'][1])),flush=True)
    synthetic=[]
    for xx in [-.1,.5,1.1]:
        for d in [1e-2,1e-6,1e-10,1e-14]:
            p=np.array([[xx,d,0.]]);a=np.array([[0.,0.,0.]]);b=np.array([[1.,0.,0.]])
            stable,_=stable_bound(p,a,b,np.array([[False]]));ref=np.array([float(v) for v in decimal_pair(p[0],a[0],b[0],bound_only=True)])
            error=float(norm(stable[0,0]-ref));threshold=1e-10+1e-10*float(norm(ref))
            native=native_parts(p,a,b,np.array([[False]]))[0][0,0]
            synthetic.append(dict(x=xx,d=d,stable=stable[0,0].tolist(),native=native.tolist(),reference=ref.tolist(),stable_error=error,native_error=float(norm(native-ref)),eligible=error<=threshold))
    limits=[]
    for x in [-.1,.5,1.1]:
        try:
            value,_=stable_bound(np.array([[x,0.,0.]]),np.array([[0.,0.,0.]]),np.array([[1.,0.,0.]]),np.array([[False]]))
            eligible=x!=.5 and np.all(value==0);status='zero'
        except ValueError:eligible=x==.5;status='contact'
        limits.append(dict(x=x,status=status,eligible=bool(eligible)))
    checks=dict(probe_unchanged=all(r['original_outputs_unchanged'] for r in parent['runs']),
        native_velocity=all(r['native_velocity_reproduction']<=1e-10 for r in RECORDS),
        native_force=all(r['native_force_reproduction']<=1e-10 for r in RECORDS),
        decimal_self=all(v['eligible'] for r in RECORDS for v in r['decimal_comparisons']),
        stable_pairs=all(r['stable_accuracy_eligible'] for r in RECORDS),
        mirror_inputs=all(r['geometry_mirror_error']<=1e-9 and r['gamma_mirror_error']<=1e-9 for r in RECORDS),
        stable_lateral=all(abs(r['stable_coefficients']['force'][1])<=1e-9 and abs(r['stable_coefficients']['moment'][2])<=1e-9 for r in RECORDS),
        bound_only_explains=all(r['bound_replacement_explains_lateral'] for r in RECORDS),
        synthetic=all(r['eligible'] for r in synthetic),limits=all(r['eligible'] for r in limits))
    result=dict(protocol='docs/verification/avl-collinear-kernel-experiment.md',source_sha256=digest(Path(__file__)),
        dependency_sha256={'check_q9_plate_convergence.py':digest(Path(memory.__file__))},parent_result_sha256=digest(parent_path),
        cases=RECORDS,synthetic=synthetic,limits=limits,checks=checks,all_criteria_pass=all(checks.values()),
        seconds=time.monotonic()-START,peak_process_working_bytes=memory.memory_bytes(),
        limits_description='postprocessing with unchanged gamma; no repaired Fortran reference or full load acceptance')
    (ROOT/'docs/verification/avl-collinear-kernel-result.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(checks=checks,seconds=result['seconds'],peak_bytes=result['peak_process_working_bytes'])),flush=True)
    if not all(checks.values()):raise SystemExit(1)

if __name__=='__main__':
    try:main()
    except Exception as exc:
        (ROOT/'build/analysis/avl-collinear-kernel-stopped.json').write_text(json.dumps(dict(error=repr(exc),cases=RECORDS),indent=2)+'\n')
        raise
