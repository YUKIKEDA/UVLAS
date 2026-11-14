"""Independent target jets, physical inertia, and restart on the same full system."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import copy
import hashlib
import json
import numpy as np
import pg_ritz_targets as target
import pg_ritz_mechanics as mechanics
from pg_ritz_state_reference import Reference, digest

ROOT=rd_paths.ROOT


def main():
    dest=ROOT/'rd/structure/main-shell/paired-geometry/records/pg-ritz-independent-state-result.json'
    assert not dest.exists()
    records=[]
    p=np.array([[.173,.251,-.72],[.647,.381,.13],[.913,.822,.88]])
    for parameters in [(.4,0.,0.),(.4,.2,0.),(.4,0.,.07),(.4,-.2,.07)]:
        for h in [.01,.0001]:
            def independent(p):
                v,d,_=target.jets(p[:,:2],parameters)
                H=h*(1+.1*p[:,0]+.05*p[:,1]);z=H*p[:,2]/2
                r,A,n,dn,g=target.geometry.graph(p[:,:2],parameters)
                return r+z[:,None]*n,v[0]+z[:,None,None]*d
            DX=[];DU=[]
            for a in range(3):
                q=p.astype(complex);q[:,a]+=1j*1e-30
                X,U=independent(q);DX.append(X.imag/1e-30);DU.append(U.imag/1e-30)
            J=np.stack(DX,axis=2)
            G=np.einsum('pati,pij->patj',np.stack(DU,axis=3),np.linalg.inv(J))
            E=(G+G.transpose(0,3,2,1))/2
            M=np.stack((E[:,0,:,0],E[:,1,:,1],E[:,2,:,2],np.sqrt(2)*E[:,1,:,2],np.sqrt(2)*E[:,0,:,2],np.sqrt(2)*E[:,0,:,1]),axis=1)
            actual=target.fields(p,h,parameters,True)[0]
            error=float(np.linalg.norm(M-actual)/np.linalg.norm(actual))
            assert error<1e-9
            records.append(dict(case='independent_target_complex_jets',parameters=parameters,h=h,error=error))
    parent=ROOT/'rd/structure/main-shell/paired-geometry/records/pg-ritz-mechanics-v2-result.json'
    report=json.loads(parent.read_text());state_file=ROOT/report['state_artifact']['artifact']
    assert hashlib.sha256(state_file.read_bytes()).hexdigest()==report['state_artifact']['artifact_sha256']
    a=np.load(state_file)
    model=mechanics.shell.prepare(1,(.4,.2,.05));data=mechanics.prepare(model)
    rigid=np.zeros((model['ns'],3,6));rigid[:model['nr'],:,:3]=np.eye(3)[None]
    physical=np.zeros((len(data['X']),3,6));physical[:,:,:3]=np.eye(3)[None]
    for j in range(3):
        rigid[:,:,3+j]=np.cross(np.eye(3)[j],a['q0'])
        physical[:,:,3+j]=np.cross(np.eye(3)[j],data['X'])
    Q=rigid.reshape(3*model['ns'],6)
    direct=np.einsum('p,pai,paj->ij',data['weights']*data['rho'],physical,physical)
    inertia=float(np.linalg.norm(Q.T@a['fullM']@Q-direct)/np.linalg.norm(direct))
    assert inertia<1e-9
    records.append(dict(case='inertia_from_analytic_physical_X',error=inertia))
    identity=hashlib.sha256(state_file.read_bytes()+parent.read_bytes()).hexdigest()
    ref=Reference(a['M'],a['K'],a['modes'],a['omega2'],identity)
    k=np.arange(1,len(a['omega2'])+1)
    amplitudes=1e-5/k**2
    initial=ref.state(a['modes']@amplitudes,a['modes']@(amplitudes*np.sqrt(a['omega2'])*.31))
    def energy(s):
        x=np.array(s['x']);v=np.array(s['v'])
        return (x@a['K']@x+v@a['M']@v)/2
    original=digest(initial);trial=ref.trial(initial,.0003)
    assert digest(initial)==original
    accepted=ref.commit(initial,trial)
    rejected=0
    for state in [initial,accepted]:
        try:ref.commit(state,trial)
        except ValueError:rejected+=1
    # A separate unconsumed candidate must also reject a different parent generation.
    stale=ref.trial(initial,.0003)
    try:ref.commit(accepted,stale)
    except ValueError:rejected+=1
    continuous=copy.deepcopy(initial)
    for _ in range(12):continuous=ref.commit(continuous,ref.trial(continuous,.0003))
    resumed=copy.deepcopy(initial)
    for _ in range(5):resumed=ref.commit(resumed,ref.trial(resumed,.0003))
    resumed=ref.load(ref.save(resumed))
    for _ in range(7):resumed=ref.commit(resumed,ref.trial(resumed,.0003))
    assert digest(continuous)==digest(resumed)
    eerror=abs(energy(continuous)-energy(initial))/energy(initial)
    assert eerror<1e-8
    for key,value in [('model','wrong'),('method','wrong'),('dtype','float32'),('dofs',0),('history_version',2),('history',{'acceleration':[]}),('accepted',False),('step',-1),('cursor',99),('x',[]),('v',[float('inf')]),('time',float('nan'))]:
        bad=copy.deepcopy(initial);bad[key]=value
        try:ref.validate(bad)
        except ValueError:rejected+=1
        else:raise AssertionError('accepted corrupt field '+key)
    wrapper=json.loads(ref.save(initial));wrapper['payload']['x'][0]+=.1
    try:ref.load(json.dumps(wrapper))
    except ValueError:rejected+=1
    else:raise AssertionError('accepted corrupt checksum')
    assert digest(initial)==original and rejected==16
    records.append(dict(case='full_DOF_exact_modal_restart_reference',dofs=len(k),all_modes_excited=True,bitwise_replay=True,energy_error=eerror,negative_controls=rejected,product_integrator_adopted=False))
    sources=[Path(__file__),ROOT/'rd/structure/main-shell/paired-geometry/experiments/pg_ritz_state_reference.py',parent,state_file,Path(target.__file__),Path(mechanics.__file__)]
    result=dict(passed=True,records=records,R4_complete=False,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    dest.write_text(json.dumps(result,indent=2));print(json.dumps(result['records'][-2:],indent=2))


if __name__=='__main__':main()
