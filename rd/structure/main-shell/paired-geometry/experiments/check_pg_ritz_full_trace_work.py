"""All outer-face and unequal-basis interface work; no multipatch PDE claim."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json
import numpy as np
import pg_ritz_shell_v3 as shell
import pg_ritz_mechanics as mechanics

ROOT=rd_paths.ROOT


def face(model,axis,position,sign,h=.01,order=17):
    n=model['n'];g,w=np.polynomial.legendre.leggauss(order)
    edge=np.concatenate([(j+(g+1)/2)/(2*n) for j in range(2*n)])
    ew=np.tile(w/(4*n),2*n)
    if axis==2:
        plane,weight=shell.lift.geometry.quadrature(n,order)
        p=np.column_stack((plane,np.full(len(plane),position)))
    else:
        z,wz=np.polynomial.legendre.leggauss(9)
        p=np.zeros((len(edge)*len(z),3));p[:,axis]=position
        p[:,1-axis]=np.repeat(edge,len(z));p[:,2]=np.tile(z,len(edge))
        weight=(ew[:,None]*wz).ravel()
    N,G,B,det,X=shell.kinematics(model,p,h,True)
    _,J,_,H,Ha=shell.lift.geometry.reference(p[:,:2],p[:,2],h,model['parameters'],True)
    area=sign*det[:,None]*np.linalg.inv(J)[:,axis,:]
    return p,weight,N,X,area


def main():
    dest=ROOT/'rd/structure/main-shell/paired-geometry/records/pg-ritz-full-trace-work-result.json';assert not dest.exists()
    rng=np.random.default_rng(20261105);models=[shell.prepare(n,(.4,-.2,.07)) for n in [1,2]]
    records=[];sigma=np.array([[.3,.04,-.02],[.04,.2,.03],[-.02,.03,.1]])
    for m in models:
        q0=mechanics.reference_coefficients(m);q=rng.normal(size=(m['ns'],3))
        force=np.zeros(3);moment=np.zeros(3);worst=0.;maperror=0.
        for axis in range(3):
            for sign in [-1,1]:
                side=sign if axis==2 else (sign+1)/2
                p,w,N,X,area=face(m,axis,side,sign)
                traction=area@sigma.T
                f=N.T@(w[:,None]*traction)
                direct=np.sum(w[:,None]*(N@q)*traction)
                error=abs(np.sum(q*f)-direct)/max(1.,abs(direct))
                worst=max(worst,float(error));maperror=max(maperror,float(np.max(abs(N@q0-X))))
                force+=np.sum(w[:,None]*traction,axis=0)
                moment+=np.sum(w[:,None]*np.cross(X,traction),axis=0)
        assert max(worst,maperror,np.linalg.norm(force),np.linalg.norm(moment))<1e-9
        records.append(dict(case='six_outer_faces',n=m['n'],work=worst,reference_map=maperror,net_force=float(np.linalg.norm(force)),net_moment=float(np.linalg.norm(moment))))
    # Common independent quadrature of a physical interface; different scalar bases.
    p,w,N2,X,area=face(models[1],0,.5,1)
    N1=shell.kinematics(models[0],p,.01,True)[0]
    q1=rng.normal(size=(models[0]['ns'],3));q2=rng.normal(size=(models[1]['ns'],3))
    lam=np.column_stack((1+p[:,1],.3-p[:,2],p[:,1]*p[:,2]))
    measure=w*np.linalg.norm(area,axis=1)
    f1=N1.T@(measure[:,None]*lam);f2=-N2.T@(measure[:,None]*lam)
    direct=np.sum(measure[:,None]*lam*(N1@q1-N2@q2))
    dual=abs(np.sum(q1*f1)+np.sum(q2*f2)-direct)/max(1.,abs(direct))
    q01=mechanics.reference_coefficients(models[0]);q02=mechanics.reference_coefficients(models[1])
    affine=np.array([[.03,.01,-.02],[.01,-.01,.02],[.01,.02,.01]])
    mismatch=float(np.max(abs(N1@(q01@affine)-N2@(q02@affine))))
    # A d-only or c-only mismatch is invisible at zeta=0 and visible on full trace.
    controls=[];m=models[1]
    for offset,count in [(m['nr'],m['nd']),(m['nr']+m['nd'],m['nc'])]:
        index=offset+np.argmax(np.linalg.norm(N2[:,offset:offset+count],axis=0))
        q=np.zeros((m['ns'],3));q[index,0]=1
        mid=p.copy();mid[:,2]=0
        invisible=float(np.max(abs(shell.kinematics(m,mid,.01,True)[0]@q)))
        detected=float(np.max(abs(N2@q)))
        assert invisible<1e-12 and detected>1e-10
        controls.append(dict(midplane=invisible,full_trace=detected))
    # Reverse the normal coordinate and the director coefficient together.
    flipped=p.copy();flipped[:,2]*=-1
    F=shell.fields(m,p[:,:2])[0]
    H=.01*(1+.1*p[:,0]+.05*p[:,1]);z=H*p[:,2]/2
    reversed_N=F[0]+(-z[:,None])*(-F[1])+z[:,None]**2*F[2]/2
    reversal=float(np.linalg.norm(reversed_N-N2)/np.linalg.norm(N2))
    assert max(dual,mismatch,reversal)<1e-9
    records.append(dict(case='unequal_basis_interface_operator',orders=[1,2],dual_work=float(dual),affine_shared_field=mismatch,normal_coordinate_reversal=reversal,midplane_negative_controls=controls,multipatch_equilibrium_verified=False))
    paths=[Path(__file__),Path(shell.__file__),Path(mechanics.__file__)]
    result=dict(passed=True,records=records,R4_complete=False,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    dest.write_text(json.dumps(result,indent=2));print(json.dumps(records,indent=2))


if __name__=='__main__':main()
