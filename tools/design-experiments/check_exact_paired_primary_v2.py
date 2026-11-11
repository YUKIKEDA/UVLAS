"""Full primary rejection/acceptance evidence for the coefficient PG candidate."""
from pathlib import Path
import gc,hashlib,json,time
import numpy as np
import exact_paired_protected as shell
import exact_paired_finite as engine
engine.shell=shell
ROOT=Path(__file__).resolve().parents[2]


def local_equivalence(model,h,variable):
    plane=np.random.default_rng(20261105).uniform(.001,.999,(19,2));z,wz=np.polynomial.legendre.leggauss(7)
    F=shell.fields(model,plane);Q,R,scales,weight,p=engine.universal(model,plane,z,wz,h,variable)
    jets=np.zeros((len(plane),24,3*model['ns']))
    values=[F[1][0],F[2][0],F[0][1],F[1][1],F[2][1],F[0][2],F[1][2],F[2][2]]
    for a in range(8):
        for j in range(3):jets[:,3*a+j,j::3]=values[a]/scales[3*a+j]
    N,G,B,det,X=shell.kinematics(model,p,h,variable)
    direct=(engine.CH@B).reshape(len(plane),len(z),6,-1)*weight[:,:,None,None]
    compressed=Q@(R@jets)
    error=float(np.linalg.norm(direct.reshape(compressed.shape)-compressed)/np.linalg.norm(direct))
    M=[]
    for order in [4,7]:
        zz,ww=np.polynomial.legendre.leggauss(order)
        pp,N,G,det,X=engine.kinematics_from_fields(model,plane,zz,h,variable,F)
        M.append(N.T@((np.tile(ww,len(plane))*det)[:,None]*N))
    scale=1/np.sqrt(np.diag(M[1]));ma=(M[0]-M[1])*scale[:,None]*scale[None,:]
    mb=M[1]*scale[:,None]*scale[None,:]
    mass=float(np.linalg.norm(ma)/np.linalg.norm(mb))
    return dict(full_3d_operator=error,mass_4_versus_7=mass)


def coarse_inclusion(model):
    n=model['n'];p,w=shell.geometry.quadrature(n,17);F=shell.fields(model,p)
    scalar=model['coarse_coordinates'];q=np.zeros((model['ns'],3,8))
    for j in range(3):q[:,j]=scalar[:,j*8:(j+1)*8]
    from pg_ritz_coarse import transform
    T=transform(n);v,d,ds=shell.compatible_particular(model['lift'],p)
    v=[a[:,:,:32*n]@T for a in v];d=d[:,:,:32*n]@T;ds=[a[:,:,:32*n]@T for a in ds]
    jets=[]
    for a in range(3):
        r=(F[a][0]@q.reshape(model['ns'],24)).reshape(len(p),3,8)
        dd=(F[a][1]@q.reshape(model['ns'],24)).reshape(len(p),3,8)
        targetd=d if a==0 else ds[a-1]
        jets.append(dict(r=float(np.linalg.norm(r-v[a])/(np.linalg.norm(v[a]) if np.linalg.norm(v[a]) else 1.)),d=float(np.linalg.norm(dd-targetd)/(np.linalg.norm(targetd) if np.linalg.norm(targetd) else 1.))))
    return jets


def main():
    dest=ROOT/'docs/verification/exact-paired-primary-v2-result.json';assert not dest.exists()
    out=ROOT/'build/analysis'/('exact-paired-primary-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir()
    parent=ROOT/'docs/verification/exact-paired-primary-result.json';previous=json.loads(parent.read_text());records=list(previous['records']);assert len(records)==4;stops=[];start=time.monotonic()
    for parameters,n,variable in [((.4,0.,0.),2,False),((.4,0.,0.),4,False),((.4,0.,.07),2,False),((.4,0.,.07),4,False),((.4,-.2,.07),2,True),((.4,-.2,.07),4,True)]:
        try:
            tick=time.monotonic();model=shell.prepare(n,parameters);prep=time.monotonic()-tick
            cold_prep=prep+(model['canonical']['seconds'] if model['canonical'].get('runtime_cache_hit',False) else 0.)
            inclusion=coarse_inclusion(model)
            if max(max(v.values()) for v in inclusion)>1e-8:raise ValueError('coarse_inclusion_jet')
            print('PREPARED',parameters,n,model['ns'],cold_prep,inclusion,flush=True)
            for h in [.01,.0001]:
                tick=time.monotonic();equivalence=local_equivalence(model,h,variable)
                if max(equivalence.values())>1e-12:raise ValueError('thickness_row_equivalence')
                sol=engine.solve(model,h,variable);solve_seconds=time.monotonic()-tick
                print('SOLVED',parameters,n,h,solve_seconds,flush=True)
                o=engine.observe(model,h,sol,variable);seconds=time.monotonic()-tick
                face=np.array([[0.,y,z] for y in np.linspace(0,1,113) for z in [-1.,-.31,0.,.47,1.]])
                N,G,B,det,X=shell.kinematics(model,face,h,variable)
                Eg,Ug=engine.coarse.fields(model,shell,face,h,variable)
                u=(N@sol['coef'].reshape(model['ns'],9)).reshape(len(face),3,3)+Ug@sol['beta']
                fixed=float(np.max(abs(u)));trace=float(np.max(abs(N)))
                checks=dict(stationarity=max(o['stationarity'])<=1e-8,rcond=sol['rcond']>=1e-10,mass=sol['mass_min']>1e-10,
                    fixed=fixed<=1e-9,trace=trace<=1e-9,resource=cold_prep+seconds<=300.)
                artifact=out/f'p{parameters}-n{n}-h{h:g}.npz'
                np.savez_compressed(artifact,coef=sol['coef'],beta=sol['beta'],scale=sol['scale'],root=sol['root'],massroot=sol['massroot'],
                    Rcoeff=model['Rcoeff'],Dcoeff=model['Dcoeff'],Tc=model['Tc'],Cindices=model['Cindices'],
                    dimensions=np.array([model['nr'],model['nd'],model['nc']]),coarse_coordinates=model['coarse_coordinates'])
                r=dict(parameters=parameters,n=n,h=h,variable=variable,ns=model['ns'],observation=o,equivalence=equivalence,
                    coarse_inclusion=inclusion,rcond=sol['rcond'],mass_min=sol['mass_min'],mass_root_ranks=sol['mass_root_ranks'],stiffness_ranks=sol['ranks'],
                    fixed=fixed,trace=trace,checks=checks,cold_preparation=cold_prep,solve_seconds=solve_seconds,
                    seconds=seconds,cold_case_seconds=cold_prep+seconds,canonical=model['canonical'],
                    artifact=artifact.relative_to(ROOT).as_posix(),artifact_sha256=hashlib.sha256(artifact.read_bytes()).hexdigest())
                records.append(r);(out/'checkpoint.json').write_text(json.dumps(dict(records=records,stops=stops),indent=2),encoding='utf-8')
                print(r,flush=True);del sol;gc.collect()
                if not all(checks.values()):raise ValueError('primary_numerical_gate')
        except Exception as exc:
            stops.append(dict(parameters=parameters,n=n,reason=repr(exc)));print('STOP',stops[-1],flush=True);break
    paths=[parent,Path(__file__),Path(shell.__file__),Path(engine.__file__),ROOT/'tools/design-experiments/exact_paired_ancf.py',ROOT/'tools/design-experiments/paired_exact_polynomials.py',ROOT/'tools/design-experiments/paired_coarse_polynomials.py',ROOT/'docs/verification/exact-paired-r4-protocol.md']
    result=dict(all_cases=len(records)==16 and not stops,records=records,stops=stops,R4_complete=False,main_shell_adopted=False,
                seconds=time.monotonic()-start,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    dest.write_text(json.dumps(result,indent=2),encoding='utf-8')


if __name__=='__main__':main()
