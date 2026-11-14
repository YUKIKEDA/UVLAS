"""Read-only two-configuration old/new kernel and first-map comparison."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,time
import numpy as np
from scipy.linalg import qr
import mq3_factor_space as old
import mq3i1_shared_space as new

ROOT=rd_paths.ROOT


def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    dest=ROOT/'rd/structure/main-shell/mq3/records/mq3i1-shared-operator-comparison.json'
    if dest.exists():raise FileExistsError(dest)
    start=time.monotonic()
    op=ROOT/'rd/structure/main-shell/mq3/records/mq3-factor-v2-result.json';npth=ROOT/'rd/structure/main-shell/mq3/records/mq3i1-shared-result.json'
    parent=json.loads(op.read_text());trial=json.loads(npth.read_text())
    assert parent['all_checks_pass'] and trial['all_checks_pass']
    out=ROOT/'build/analysis'/('mq3i1-compare-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir()
    (out/'old').mkdir();(out/'new').mkdir()
    old.family.configure('covariant',True);new.family.configure('covariant',True)
    caches=[]
    for source,family in [(parent,old.family),(trial,new.family)]:
        cache={}
        keys=['q','T','rows','U'] if family is old.family else new.family.KEYS
        for key,item in source['maps'].items():
            assert digest(ROOT/item['artifact'])==item['artifact_sha256']
            saved=np.load(ROOT/item['artifact'])
            cache[key]=dict(state={k:saved[k] for k in keys},X0=saved['X0'],Xd=saved['Xd'],
                validation=item['validation'],metrics=item['metrics'],existing=True)
        caches.append(cache)
    records=[];stop=None
    try:
        for amp in [0.,.2]:
            ro=next(r for r in parent['records'] if (r['metric'],r['work'],r['kappa'],r['amplitude'],r['n'])==('covariant',True,.4,amp,2))
            rn=next(r for r in trial['records'] if (r['kappa'],r['amplitude'])==(.4,amp))
            for r in [ro,rn]:assert digest(ROOT/r['artifact'])==r['artifact_sha256']
            so,sn=np.load(ROOT/ro['artifact']),np.load(ROOT/rn['artifact'])
            ao=old.assemble(.4,2,amp,11,caches[0],{},out/'old')
            an=new.assemble(.4,2,amp,11,caches[1],{},out/'new')
            Qo,_=qr(so['kernel'],mode='economic',check_finite=False)
            Qn,_=qr(sn['kernel'],mode='economic',check_finite=False)
            projector=float(np.linalg.norm(Qo@Qo.T-Qn@Qn.T))
            cross=lambda A,N:float(np.linalg.norm(A@N)/(np.linalg.norm(A)*np.linalg.norm(N)))
            old_on_new=cross(ao['A0'],Qn);new_on_old=cross(an['A0'],Qo)
            Lo=np.column_stack((ao['A1']@Qn,ao['D1']));Ln=np.column_stack((an['A1']@Qn,an['D1']))
            first=float(np.linalg.norm(Lo-Ln)/max(1.,np.linalg.norm(Lo)))
            observed={}
            for source,saved in [('old',so),('new',sn)]:
                for kind in ['quadratic','quartic']:
                    u=saved['u_'+kind];rd=u[an['free'],:6].reshape(-1);c=u[an['free'],6:].reshape(-1)
                    eold=float(np.linalg.norm(ao['A1']@rd+ao['D1']@c-ao['targets'][kind])/np.linalg.norm(ao['targets'][kind]))
                    enew=float(np.linalg.norm(an['A1']@rd+an['D1']@c-an['targets'][kind])/np.linalg.norm(an['targets'][kind]))
                    observed[source+'_'+kind]=dict(old_error=eold,new_error=enew,difference=abs(eold-enew))
            artifact=out/f'amp{amp:g}.npz'
            np.savez_compressed(artifact,A0old=ao['A0'],A0new=an['A0'],A1old=ao['A1'],A1new=an['A1'],D1old=ao['D1'],D1new=an['D1'],Qold=Qo,Qnew=Qn)
            records.append(dict(amplitude=amp,old_nullity=Qo.shape[1],new_nullity=Qn.shape[1],
                kernel_projector_difference=projector,old_on_new_kernel=old_on_new,new_on_old_kernel=new_on_old,
                numerically_same_kernel=Qo.shape[1]==Qn.shape[1] and projector<=1e-7 and max(old_on_new,new_on_old)<=1e-9,
                first_operator_difference=first,first_operator_agrees=first<=1e-9,saved_solution_observations=observed,
                artifact=artifact.relative_to(ROOT).as_posix(),artifact_sha256=digest(artifact)))
            print(amp,records[-1],flush=True)
    except Exception as exc:stop=repr(exc)
    sources=[Path(__file__),Path(old.__file__),Path(new.__file__),op,npth,ROOT/'rd/structure/main-shell/mq3/records/mq3i1-shared-followup-protocol.md']
    dest.write_text(json.dumps(dict(records=records,complete=len(records)==2 and stop is None,stop_reason=stop,
        sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start),indent=2),encoding='utf-8')


if __name__=='__main__':main()
