"""Read-only residual diagnosis after the first shared gate failure; no new solve."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,math
import numpy as np
import reference_work_fit as fit
f=fit.f;c=f.c;ROOT=rd_paths.ROOT
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def tensor_energy(E):return .3/(1.3*.4)*np.trace(E,axis1=1,axis2=2)**2+2/2.6*np.einsum('pij,pij->p',E,E)
def analytical_flat(q,h):
    # Exact polynomial admissible quartic field on this rectangular affine cell.
    x=q.reshape(75,3)[:25,0];a=(x.max()+x.min())/2;b=(x.max()-x.min())/2;coeff=np.zeros((len(c.IDS),3))
    for n,component,factor,z in [(4,2,1/12,0),(3,0,-h/6,1)]:
        for power in range(n+1):
            mon=math.comb(n,power)*a**(n-power)*b**power
            for degree,value in c.powers(power):key=(degree,0,z);coeff[c.LOOK[key],component]+=factor*mon*float(value)/c.NORM[c.LOOK[key]]
    return coeff
def main():
    dest=ROOT/'rd/structure/main-shell/reference-work/records/reference-work-v3-residual-diagnostic.json';assert not dest.exists();parent=ROOT/'rd/structure/main-shell/reference-work/records/reference-work-v3-shared-result.json';source=json.loads(parent.read_text());records=[]
    for r in source['records']:
        assert r['kappa']==0. and r['amplitude']==0.
        data=np.load(ROOT/r['artifact']);h=r['h'];sums=np.zeros(8);reference=0.;maxcondition=0.
        for ids,values in zip(data['cells'],data['local']):
            q=data['ref'][ids].reshape(-1);s=f.prepare(q,h);u=values[:,0];star=fit.target(q,s['p'],h,0.,'quartic');norm=h**1.5
            target=(star@f.m.CH.T*np.sqrt(s['weight'])[:,None]).reshape(-1)/norm;b=np.zeros(len(s['root']));b[:len(target)]=target
            W=s['root']/norm;res=W@u-b;res128=np.asarray(np.asarray(W,dtype=np.longdouble)@np.asarray(u,dtype=np.longdouble)-np.asarray(b,dtype=np.longdouble),dtype=float)
            physical=s['B']@u-star;weighted=(physical@f.m.CH.T*np.sqrt(s['weight'])[:,None]).reshape(-1)/norm
            sums[0]+=res@res;sums[1]+=res128@res128;sums[2]+=weighted@weighted
            sums[2]+=sum(s['tau']*face['measure']@np.sum((face['jump']@u.reshape(-1,3))**2,axis=1)/h**3 for face in s['faces'])
            scale=1/np.linalg.norm(W,axis=0);R=f.qrtools.rfactor(np.column_stack((W*scale,b)));nr=W.shape[1];rr=R[:,:nr]@(u/scale)-R[:,nr];sums[3]+=rr@rr
            maxcondition=max(maxcondition,1/f.qrtools.rcond(R[:3*s['ni'],:3*s['ni']]))
            coeff=analytical_flat(q,h);oracle=np.zeros((s['nscalar'],3));oracle[:75]=coeff[f.VI]
            for face in s['faces']:oracle[face['offset']:face['offset']+face['n']]=face['Mcoeff'].T@c.restrict(coeff,face['axis'],face['sign'])
            Eoracle=s['B']@oracle.reshape(-1)-star;sums[4]+=s['weight']@np.einsum('pk,kl,pl->p',Eoracle,f.m.C,Eoracle)/h**3
            for face in s['faces']:sums[4]+=s['tau']*face['measure']@np.sum((face['jump']@oracle)**2,axis=1)/h**3
            _,g=c.scalar_basis(s['p']);du=np.einsum('pjn,na->paj',g,coeff);grad=du@s['invJ'];E=(grad+grad.transpose(0,2,1))/2;target_tensor=np.einsum('kij,pk->pij',f.m.MB,star)
            sums[5]+=s['weight']@tensor_energy(E-target_tensor)/h**3
            si=f.prepare(q,h,order=15,zorder=9,independent=True);stari=fit.target(q,si['p'],h,0.,'quartic');ei=si['B']@u-stari;sums[6]+=si['weight']@np.einsum('pk,kl,pl->p',ei,f.m.C,ei)/h**3
            for face in si['faces']:sums[6]+=si['tau']*face['measure']@np.sum((face['jump']@u.reshape(-1,3))**2,axis=1)/h**3
            # Reassociate the same unmodified weak-gradient operator.
            weak=s['Sbody']@(s['D']@u.reshape(-1,3));gradw=np.einsum('pja,pjk->pak',weak,s['invJ']);Ew=(gradw+gradw.transpose(0,2,1))/2
            sums[7]+=s['weight']@tensor_energy(Ew-target_tensor)/h**3
            for face in s['faces']:sums[7]+=s['tau']*face['measure']@np.sum((face['jump']@u.reshape(-1,3))**2,axis=1)/h**3
            reference+=target@target
        error=np.sqrt(sums/reference);names=['direct_root_residual','extended_precision_root_product','primary_tensor','rebuilt_QR_residual','exact_admissible_polynomial_weak_operator','exact_polynomial_strong_gradient','independent_tensor','reassociated_weak_gradient']
        records.append(dict(h=h,errors=dict(zip(names,map(float,error))),reported_QR_error=r['fit']['error'][0],reported_independent_error=r['observed']['error'][0],local_condition_estimate=maxcondition))
    ps=[Path(__file__),Path(fit.__file__),ROOT/'rd/structure/main-shell/reference-work/experiments/reference_work_ancf_v3.py',parent,ROOT/'rd/structure/main-shell/reference-work/records/reference-work-v3-protocol.md']
    dest.write_text(json.dumps(dict(records=records,new_global_solve=False,acceptance_thresholds_changed=False,failed_gate_overridden=False,unique_cause_proven=False,sources={p.relative_to(ROOT).as_posix():digest(p) for p in ps}),indent=2),encoding='utf-8');print(json.dumps(records,indent=2))
if __name__=='__main__':main()
