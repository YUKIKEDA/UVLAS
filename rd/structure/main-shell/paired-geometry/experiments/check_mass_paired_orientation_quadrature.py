"""Independent quadrature, all face orientations, and layer reversal."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,time,traceback,gc
import numpy as np
import mass_paired_mechanics as m
ROOT=rd_paths.ROOT
def rel(a,b):return float(np.linalg.norm(a)/max(1e-30,np.linalg.norm(b)))

def main():
    dest=ROOT/'rd/structure/main-shell/paired-geometry/records/mass-paired-orientation-result.json';assert not dest.exists();records=[];stops=[];start=time.monotonic()
    try:
        model=m.shell.prepare(1,(.4,.2,.07,.1,-.06,.025),.05,True,'extended');q0=model['reference_coordinates'];rng=np.random.default_rng(20261105)
        direction=rng.normal(size=q0.shape);sample=np.random.default_rng(16).uniform([0.,0.,-1.],[1.,1.,1.],(257,3))
        G=m.shell.kinematics(model,sample,.05,True)[1];direction*=.003/np.max(abs(np.einsum('na,pnb->pab',direction,G)))
        u=q0@(m.FSTAR-np.eye(3)).T+direction
        previous=None
        for order in [45,49]:
            data=m.prepare(model,.05,order=order,zorder=7);U,f,K,_=m.evaluate(data,u,True)
            current=dict(U=U,f=f,K=K,M=data['M'])
            if previous is not None:
                metrics={key:rel(np.asarray(current[key])-np.asarray(previous[key]),np.asarray(current[key])) for key in current}
                assert max(metrics.values())<=1e-8;records.append(dict(case='independent_nonlinear_volume_quadrature',orders=[45,49],metrics=metrics))
            previous=current;del data;gc.collect()
        shell=m.shell;h=.05;plane,wp=shell.geometry.quadrature(1,25);q,w=np.polynomial.legendre.leggauss(7)
        masses=[]
        # Reverse thickness orientation and the associated layer/density ordering.
        for reverse,flip_layers in [(False,False),(True,True),(True,False)]:
            M=np.zeros((model['ns'],model['ns']))
            for lo,hi,rho in [(-1.,0.,1.2),(0.,1.,.8)]:
                zz=(lo+hi)/2+(hi-lo)*q/2;p=np.column_stack((np.repeat(plane,len(q),axis=0),np.tile(zz,len(plane))))
                old=p.copy()
                if reverse:old[:,2]*=-1
                N,G,B,det,X=shell.kinematics(model,old,h,True)
                density=2-rho if flip_layers else rho
                ww=(wp[:,None]*w*(hi-lo)/2).ravel()*det*density
                M+=N.T@(ww[:,None]*N)
            masses.append(M)
        layer=rel(masses[1]-masses[0],masses[0]);negative=rel(masses[2]-masses[0],masses[0])
        assert layer<=1e-10 and negative>1e-5
        records.append(dict(case='reversed_normal_layer_order',mass_error=layer,forgotten_layer_flip_defect=negative))
        qcart=rng.normal(size=(model['ns'],3));errors=[];negative_errors=[]
        for axis in range(3):
            for sign in [-1,1]:
                pp=np.random.default_rng(93+2*axis+sign).uniform([0.,0.,-1.],[1.,1.,1.],(41,3));pp[:,axis]=sign if axis==2 else (sign+1)/2
                old=pp.copy();old[:,2]*=-1
                X,J,det,H,Ha=shell.geometry.reference(old[:,:2],old[:,2],h,model['parameters'],True)
                Jreverse=J.copy();Jreverse[:,:,2]*=-1
                area=sign*abs(np.linalg.det(Jreverse))[:,None]*np.linalg.inv(Jreverse)[:,axis,:]
                oldsign=-sign if axis==2 else sign
                expected=oldsign*det[:,None]*np.linalg.inv(J)[:,axis,:]
                P=np.zeros((len(X),3,3));P[:,0,0]=1+X[:,1]**2;P[:,1,2]=X[:,0]-.3;P[:,2,1]=.7+X[:,2]
                N=shell.kinematics(model,old,h,True)[0];virtual=N@qcart
                work=np.einsum('pa,pab,pb->p',virtual,P,area);reference=np.einsum('pa,pab,pb->p',virtual,P,expected)
                errors.append(rel(work-reference,reference));negative_errors.append(rel(-work-reference,reference))
        assert max(errors)<=1e-10 and min(negative_errors)>1
        records.append(dict(case='all_six_faces_reversed_chart_virtual_work',errors=errors,wrong_signed_measure_defects=negative_errors))
    except Exception as exc:stops.append(dict(reason=repr(exc),traceback=traceback.format_exc()));print('FAIL',stops[-1],flush=True)
    paths=[Path(__file__),Path(m.__file__),Path(m.shell.__file__),Path(m.kernel.__file__)]
    result=dict(passed=len(records)==3 and not stops,records=records,stops=stops,seconds=time.monotonic()-start,R4_complete=False,
        sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print('ORIENTATION_RESULT',result['passed'],flush=True)

if __name__=='__main__':main()
