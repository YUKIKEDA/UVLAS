"""Independent jets, true thickness, full trace and global reference bounds."""
from pathlib import Path
import hashlib,json,time,traceback
import numpy as np
import mass_paired_ancf as shell
import mass_paired_mechanics as mechanics
import general_paired_fields as target
from paired_reference_admissibility import certify
ROOT=Path(__file__).resolve().parents[2]
def rel(a,b):return float(np.linalg.norm(a)/max(1.,np.linalg.norm(b)))

def main():
    dest=ROOT/'docs/verification/mass-paired-reference-result.json';assert not dest.exists()
    records=[];stops=[];start=time.monotonic()
    try:
        p=np.array([[.173,.251,-.72],[.647,.381,.13],[.913,.822,.88]])
        for pars in [(.4,0.,0.),(.4,.2,0.),(.4,0.,.07),(.4,-.2,.07),(.4,.2,.07,.1,-.06,.025),(.4,.2,.07,.1,-.06,.025,0.,1.,.2)]:
            for h in [.05,.01,.0001]:
                proof=certify(pars,h,True);X,J,det,H,Ha=shell.geometry.reference(p[:,:2],p[:,2],h,pars,True)
                DX=[];DU=[];Ddv=[];v,d,ds=target.jets(p[:,:2],pars)
                for a in range(3):
                    q=p.astype(complex);q[:,a]+=1j*1e-30
                    Xc=shell.geometry.reference(q[:,:2],q[:,2],h,pars,True)[0]
                    vc,dc,_=target.jets(q[:,:2],pars);Hc=h*(1+.1*q[:,0]+.05*q[:,1]);Uc=vc[0]+(Hc*q[:,2]/2)[:,None,None]*dc
                    DX.append(Xc.imag/1e-30);DU.append(Uc.imag/1e-30)
                    if a<2:Ddv.append([vc[a+1].imag/1e-30,dc.imag/1e-30])
                JJ=np.stack(DX,axis=2);G=np.einsum('pati,pij->patj',np.stack(DU,axis=3),np.linalg.inv(JJ))
                E=(G+G.transpose(0,3,2,1))/2
                M=np.stack((E[:,0,:,0],E[:,1,:,1],E[:,2,:,2],np.sqrt(2)*E[:,1,:,2],np.sqrt(2)*E[:,0,:,2],np.sqrt(2)*E[:,0,:,1]),axis=1)
                actual=target.fields(p,h,pars,True)[0]
                metrics=dict(reference_jacobian=rel(J-JJ,J),target_physical_strain=rel(M-actual,actual),
                    v_ss=rel(Ddv[0][0]-v[3],v[3]),v_tt=rel(Ddv[1][0]-v[5],v[5]),
                    d_s=rel(Ddv[0][1]-ds[0],ds[0]),d_t=rel(Ddv[1][1]-ds[1],ds[1]))
                r,A,n,dn,g=shell.geometry.graph(p[:,:2],pars)
                metrics['unit_normal']=float(np.max(abs(np.sum(n*n,axis=1)-1)))
                metrics['real_normal_thickness']=rel(X-r-H[:,None]*p[:,2,None]*n/2,X)
                eta_ss=np.einsum('pa,pat->pt',A[:,:,0],v[1]);eta_st=(np.einsum('pa,pat->pt',A[:,:,0],v[2])+np.einsum('pa,pat->pt',A[:,:,1],v[1]))/2
                metrics['membrane_ss_st']=max(float(np.max(abs(eta_ss))),float(np.max(abs(eta_st))))
                metrics['transverse_constraints']=max(float(np.max(abs(np.einsum('pa,pat->pt',n,d)))),
                    *[float(np.max(abs(np.einsum('pa,pat->pt',n,v[a+1])+np.einsum('pa,pat->pt',A[:,:,a],d)))) for a in range(2)])
                assert max(metrics.values())<1e-8 and np.min(det)>=proof['global_det_J_lower']
                records.append(dict(case='analytic_jets_and_reference',parameters=pars,h=h,metrics=metrics,certificate=proof))
        model=shell.prepare(1,(.4,.2,.07,.1,-.06,.025),.05,True,'extended')
        N,G,B,det,X=shell.kinematics(model,p,.05,True);q0=model['reference_coordinates']
        metrics=dict(reference_position=rel(N@q0-X,X),reference_gradient=rel(np.einsum('na,pnb->pab',q0,G)-np.eye(3),np.eye(3)))
        DN=[]
        for a in range(3):
            pp=p.astype(complex);pp[:,a]+=1j*1e-30
            DN.append(shell.kinematics(model,pp,.05,True)[0].imag/1e-30)
        J=shell.geometry.reference(p[:,:2],p[:,2],.05,model['parameters'],True)[1]
        metrics['physical_basis_derivatives']=rel(np.stack(DN,axis=2)@np.linalg.inv(J)-G,G)
        Z,ranks,defect=mechanics.fixed_trace_kernel(model);metrics['trace_holdout']=defect
        assert max(metrics.values())<1e-8
        # Exact kernel identity V=S+LTS, H=(I-LT)S: compare exact coefficient ranks.
        rows,Q,_,full=shell.source_basis(1,model['parameters'],'extended')
        rows,HQ,_,hom=shell.source_basis(1,model['parameters'],'clamped')
        full_dimension=Q.ncols()+4;hom_dimension=HQ.ncols()+2
        assert full_dimension-hom_dimension==ranks[0] and Z.shape[1]==hom_dimension
        records.append(dict(case='same_full_basis_reference_and_trace',metrics=metrics,full_dimension=full_dimension,
            homogeneous_dimension=hom_dimension,trace_ranks=ranks,exact_kernel_dimension_identity=True))
        # Global Cartesian frame rotates q and the reference covectors together.
        angle=.63;R=np.array([[np.cos(angle),-np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1.]])
        rng=np.random.default_rng(18);q=rng.normal(size=(model['ns'],3));Gc=G@R.T
        original=np.einsum('na,pnb->pab',q,G);rotated=np.einsum('na,pnb->pab',q@R.T,Gc)
        frame=rel(rotated-R@original@R.T,original);assert frame<1e-12
        records.append(dict(case='fixed_global_reference_frame_covariance',error=frame))
        rejected=0
        for pars,h in [(model['parameters'],-.01),(model['parameters'],0.),(model['parameters'],100.),((.4,.2,0.,0.,0.,0.,0.,-1.),.01)]:
            try:certify(pars,h,True)
            except ValueError:rejected+=1
        assert rejected==4
        records.append(dict(case='invalid_reference_negative_controls',rejected=rejected))
    except Exception as exc:stops.append(dict(reason=repr(exc),traceback=traceback.format_exc()));print(stops[-1],flush=True)
    paths=[Path(__file__),ROOT/'tools/design-experiments/paired_reference_admissibility.py',Path(shell.__file__),Path(target.__file__),Path(shell.geometry.__file__),Path(mechanics.__file__)]
    result=dict(passed=len(records)==21 and not stops,records=records,stops=stops,seconds=time.monotonic()-start,R4_complete=False,
        sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print('REFERENCE_RESULT',result['passed'],flush=True)

if __name__=='__main__':main()
