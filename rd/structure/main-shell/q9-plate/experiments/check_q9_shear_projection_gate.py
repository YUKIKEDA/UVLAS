"""Local mixed projection candidate, gated before any global bending solve."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,itertools,json,time
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import cho_factor,cho_solve
import check_q9_distorted_plate as base
import check_q9_curved_load as geometry
import check_q9_independent_load as material
import check_q9_plate_convergence as memory
ROOT=rd_paths.ROOT;F=rd_paths.RECORDS
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def factor(M):
    d=1/np.sqrt(np.diag(M));a=d[:,None]*M*d[None,:];sv=np.linalg.svd(a,compute_uv=False)
    if np.min(sv)<=1e-12*np.max(sv):raise ValueError('rank_deficient_projection_Gram')
    f=cho_factor(a,lower=True)
    return d,a,f
def solve(f,rhs):
    d,_,c=f
    return d[:,None]*cho_solve(c,d[:,None]*rhs) if rhs.ndim==2 else d*cho_solve(c,d*rhs)

def matrices(q,h,order):
    scalar,g,weights=base.table(order,h);qg=q.reshape(27,3);jac=np.einsum('ia,pib->pab',qg,g);measure=weights*np.linalg.det(jac)
    grad=np.linalg.solve(jac.transpose(0,2,1),g.transpose(0,2,1)).transpose(0,2,1)
    B=np.zeros((len(g),6,81));B[:,0,0::3]=grad[:,:,0];B[:,1,1::3]=grad[:,:,1];B[:,2,2::3]=grad[:,:,2]
    B[:,3,1::3]=grad[:,:,2];B[:,3,2::3]=grad[:,:,1];B[:,4,0::3]=grad[:,:,2];B[:,4,2::3]=grad[:,:,0];B[:,5,0::3]=grad[:,:,1];B[:,5,1::3]=grad[:,:,0]
    normal=np.linalg.cholesky(material.elastic[:3,:3]).T
    rn=(np.sqrt(measure)[:,None,None]*np.einsum('ab,pbc->pac',normal,B[:,:3,:])).reshape(-1,81)
    rxy=np.sqrt(material.mu*measure)[:,None]*B[:,5,:];Kn=rn.T@rn+rxy.T@rxy
    Bs=B[:,[4,3],:];rs=(np.sqrt(material.mu*measure)[:,None,None]*Bs).reshape(-1,81);Ks=rs.T@rs
    p,_=leggauss(order);z,_=leggauss(3);points=np.array(list(itertools.product(p,p,z)))
    def P(t,k):return np.ones_like(t) if k==0 else (t if k==1 else (3*t*t-1)/2)
    S=np.zeros((len(points),2,36));col=0
    for component,first,second in [(0,2,3),(1,3,2)]:
        for i,j,k in itertools.product(range(first),range(second),range(3)):
            S[:,component,col]=P(points[:,0],i)*P(points[:,1],j)*P(points[:,2],k);col+=1
    V=np.linalg.solve(jac[:,:2,:2].transpose(0,2,1),S)
    M=np.einsum('p,pai,paj->ij',material.mu*measure,V,V);H=np.einsum('p,pai,paj->ij',material.mu*measure,V,Bs)
    f=factor(M);Kp=Kn+H.T@solve(f,H)
    return dict(M=M,H=H,Kn=Kn,Kfull=Kn+Ks,Kp=Kp,factor=f)

def patch_fields(q):
    nodes=q.reshape(9,3,3);fields=[]
    for component in range(6):
        E=np.zeros((3,3))
        if component<3:E[component,component]=1
        else:
            a,b=[(1,2),(0,2),(0,1)][component-3];E[a,b]=E[b,a]=.5
        fields.append((nodes@E.T).reshape(81))
    fields.append(geometry.raw_exact_q(q,False,np.eye(3),np.zeros(3)))
    return fields

def main():
    target=F/'q9-shear-projection-gate-result.json'
    if target.exists():raise FileExistsError(target)
    out=ROOT/'build/analysis'/('q9-shear-projection-gate-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir()
    start=time.monotonic();records=[];stop=None
    try:
        for name in ['RECT','SKEW','TRAP']:
            q=geometry.geometry(name);data={order:matrices(q,.1,order) for order in [7,9]};d=data[9];scale=1/np.sqrt(d['Kfull'].diagonal())
            scaled=lambda k:scale[:,None]*k*scale[None,:]
            K=scaled(d['Kp']);Kf=scaled(d['Kfull']);norm=float(np.linalg.norm(K,2));singular=np.linalg.svd(K,compute_uv=False)
            ranks={str(t):int(np.count_nonzero(singular>t*singular[0])) for t in [1e-11,1e-10,1e-9]}
            eigen=np.linalg.eigvalsh((K+K.T)/2);difference=Kf-K;difference_eigen=np.linalg.eigvalsh((difference+difference.T)/2)
            symmetry=float(np.linalg.norm(K-K.T)/np.linalg.norm(K));integration=float(np.linalg.norm(scaled(data[7]['Kp']-d['Kp']))/np.linalg.norm(Kf))
            fields=patch_fields(q);patches=[]
            for u in fields:
                denominator=np.linalg.norm(scale*(d['Kfull']@u));patches.append(float(np.linalg.norm(scale*((d['Kp']-d['Kfull'])@u))/denominator))
            rigid=[];coordinates=q.reshape(9,3,3)
            for axis in np.eye(3):
                u=np.zeros_like(coordinates);u[:,0,:]=axis;rigid.append(u.reshape(81))
                rigid.append(np.cross(np.broadcast_to(axis,coordinates.shape),coordinates).reshape(81))
            rigid_error=max(float(np.linalg.norm(K@(u/scale))/(norm*np.linalg.norm(u/scale))) for u in rigid)
            rng=np.random.default_rng(831);u=scale*rng.normal(size=81);rhs=d['H']@u;beta=solve(d['factor'],rhs);residual=float(np.linalg.norm(d['M']@beta-rhs)/np.linalg.norm(rhs))
            mixed=float(u@d['Kn']@u/2+beta@rhs-beta@d['M']@beta/2);condensed=float(u@d['Kp']@u/2);energy_error=abs(mixed-condensed)/max(1.,abs(condensed))
            sf=d['factor'][0];direct=sf[:,None]*np.linalg.solve(d['factor'][1],sf[:,None]*d['H']);Kdirect=d['Kn']+d['H'].T@direct
            condensation=float(np.linalg.norm(scaled(Kdirect-d['Kp']))/np.linalg.norm(K))
            duplicate=d['M'].copy();duplicate[-1,:]=duplicate[0,:];duplicate[:,-1]=duplicate[:,0]
            rejected=False
            try:factor(duplicate)
            except ValueError:rejected=True
            condition=float(np.linalg.cond(d['factor'][1]));checks=dict(gram=condition<=1e8 and residual<=1e-11,symmetry=symmetry<=1e-12,integration=integration<=1e-8,
                rank=all(v==75 for v in ranks.values()),rigid=rigid_error<=1e-11,positive=float(eigen[0])>=-1e-10*max(1.,norm),
                energy_order=float(difference_eigen[0])>=-1e-10*max(1.,float(np.linalg.norm(difference,2))),patch=max(patches)<=1e-10,
                mixed_energy=energy_error<=1e-10,condensation=condensation<=1e-10,duplicate_rejected=rejected)
            file=out/(name+'.npz');np.savez_compressed(file,**{k:v for k,v in d.items() if k!='factor'},q=q,rigid=np.array(rigid),patch_fields=np.array(fields))
            records.append(dict(case=name,rank=ranks,scaled_singular_values=singular.tolist(),scaled_eigenvalues=eigen.tolist(),gram_condition=condition,gram_residual=residual,
                symmetry=symmetry,integration_difference=integration,patch_force_differences=patches,rigid_action=rigid_error,minimum_energy_difference_eigenvalue=float(difference_eigen[0]),
                mixed_energy_difference=energy_error,condensation_difference=condensation,checks=checks,artifact=file.relative_to(ROOT).as_posix(),artifact_sha256=digest(file)))
            print(json.dumps(dict(case=name,ranks=ranks,checks=checks,patches=patches)),flush=True)
            if time.monotonic()-start>120 or memory.memory_bytes()>512*1024**2:raise RuntimeError('resource_limit')
        if sum(p.stat().st_size for p in out.iterdir())>16*1024**2:raise RuntimeError('output_limit')
    except Exception as exc:stop=repr(exc)
    result=dict(source_sha256=digest(__file__),candidate_sha256=digest(ROOT/'rd/structure/main-shell/q9-plate/q9-shear-projection-candidate.md'),protocol_sha256=digest(F/'q9-shear-projection-gate-experiment.md'),
        dependency_sha256={Path(m.__file__).name:digest(m.__file__) for m in [base,geometry,material,memory]},cases=records,stop_reason=stop,seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),
        all_criteria_pass=bool(stop is None and len(records)==3 and all(all(r['checks'].values()) for r in records)),product_shell_adopted=False)
    target.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if not result['all_criteria_pass']:raise SystemExit(1)
if __name__=='__main__':main()
