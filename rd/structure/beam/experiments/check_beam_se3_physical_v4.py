"""Independent rod, Timoshenko static/modal and Cosserat BVP checks."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,time,traceback
import numpy as np
from flint import arb_mat,ctx
from scipy.linalg import qr,solve_triangular,svd,cholesky,expm
from scipy.integrate import solve_bvp
from scipy.optimize import brentq
from scipy.sparse import diags,kron,eye
from scipy.sparse.linalg import spsolve
import beam_se3_v2 as b
import beam_se3_high_precision as hp
from beam_decimal_reference import sincos
from decimal import Decimal as Dec
ROOT=rd_paths.ROOT

def torsion_rectangle(width,height,terms=201):
    if height>width:width,height=height,width
    k=np.arange(1,2*terms,2,dtype=float)
    return width*height**3/3*(1-192*height/(np.pi**5*width)*np.sum(np.tanh(k*np.pi*width/(2*height))/k**5))
def section(h):
    E=1e6;nu=.3;rho=2.3;width=.2;A=width*h;Iy=A*width**2/12;Iz=A*h*h/12;G=E/(2*(1+nu));J=torsion_rectangle(width,h)
    K=np.diag([E*A,5/6*G*A,5/6*G*A,G*J,E*Iy,E*Iz]);mass=b.spatial_mass(rho*A,np.zeros(3),np.diag([rho*(Iy+Iz),rho*Iy,rho*Iz]))
    return K,mass,dict(E=E,rho=rho,A=A,Iy=Iy,Iz=Iz,J=J,G=G)
def linear_operators(n,h,needmass=False):
    L=1/n;K,mass,props=section(h);r=np.array([[0.,0.,0.],[L,0,0]]);R=np.array([np.eye(3)]*2);e=b.evaluate(r,R,L,K,np.array([1.,0,0,0,0,0]))
    local=np.linalg.cholesky(K).T@e['B']*np.sqrt(L);root=np.zeros((6*n,6*(n+1)))
    M=np.zeros((6*(n+1),6*(n+1))) if needmass else None
    if needmass:Me,_=b.inertia(r,R,L,mass,4)
    for i in range(n):
        root[6*i:6*i+6,6*i:6*i+12]=local
        if needmass:M[6*i:6*i+12,6*i:6*i+12]+=Me
    return root[:,6:],M[6:,6:] if needmass else None,props
def linear_solution(root,F):
    scales=1/np.linalg.norm(root,axis=0);Q,T=qr(root*scales,mode='economic');ss=svd(T,compute_uv=False)
    x=solve_triangular(T,solve_triangular(T.T,scales*F,lower=True));u=scales*x
    res=np.linalg.norm(T.T@(T@(u/scales))-scales*F)/np.linalg.norm(scales*F)
    if res>1e-8 or ss[-1]/ss[0]<1e-5:
        # Solve the SAME square strain root, retaining the displacement until
        # the original force residual is evaluated at high precision.
        ctx.prec=256;Ar=arb_mat(root.tolist());Fr=arb_mat([[float(x)] for x in F])
        stress=Ar.transpose().solve(Fr);state=Ar.solve(stress)
        residual=Ar.transpose()*(Ar*state)-Fr
        res=max(float(abs(residual[i,0]).upper()) for i in range(len(F)))/max(abs(F))
        u=np.array([float(state[i,0]) for i in range(len(F))])
        linear_solution.last_state=dict(precision_bits=256,coordinates=[str(state[i,0]) for i in range(len(F))],retained_residual=res)
    else:linear_solution.last_state=dict(precision_bits=53,retained_residual=float(res))
    return u,float(res),float(ss[-1]/ss[0])
def frequencies(root,M,columns):
    Mr=M[np.ix_(columns,columns)];L=cholesky(Mr,lower=True);A=solve_triangular(L,root[:,columns].T,lower=True).T
    U,s,Vt=svd(A,full_matrices=False);V=solve_triangular(L.T,Vt.T,lower=False)[:,::-1];omega=s[::-1]
    backward=np.linalg.norm(A.T@(A@Vt.T)-Vt.T*s[None,:]**2)/np.linalg.norm(A.T@(A@Vt.T))
    orth=np.linalg.norm(V.T@Mr@V-np.eye(len(columns)))/np.sqrt(len(columns))
    return omega,float(backward),float(orth)
def continuous_bending(props,h):
    K,mass,_=section(h);alpha=K[5,5]/K[1,1];beta=mass[5,5]/mass[0,0]
    def determinant(lam):
        A=np.array([[0,1,alpha,0],[0,0,0,1],[-lam*lam,0,0,0],[0,-beta*lam*lam,-1,0]],float)
        return np.linalg.det(expm(A)[2:,2:])
    roots=[];grid=np.linspace(.05,35.,1000);prev=grid[0];fp=determinant(prev)
    for x in grid[1:]:
        fx=determinant(x)
        if fx*fp<0:roots.append(brentq(determinant,prev,x,xtol=1e-12))
        if len(roots)==2:break
        prev=x;fp=fx
    if len(roots)!=2:raise ValueError('continuous modal roots')
    return np.array(roots)*np.sqrt(K[5,5]/mass[0,0])
def nonlinear(n,P,reference_K):
    L=1/n;r=np.column_stack((np.linspace(0,1,n+1),np.zeros((n+1,2))));R=np.array([np.eye(3)]*(n+1));columns=np.array([6*i+j for i in range(1,n+1) for j in [0,1,5]])
    history=[]
    for step in np.linspace(.25,1.,4):
        for iteration in range(15):
            f=np.zeros(6*(n+1));K=np.zeros((len(f),len(f)));U=0.
            for i in range(n):
                e=b.evaluate(r[i:i+2],R[i:i+2],L,reference_K,np.array([1.,0,0,0,0,0]));sl=slice(6*i,6*i+12);f[sl]+=e['f'];K[sl,sl]+=e['K'];U+=e['U']
            external=np.zeros(len(f));external[-5]=step*P;res=(f-external)[columns]
            scaling=1/np.sqrt(np.maximum(abs(np.diag(K)[columns]),1e-300));error=float(np.linalg.norm(scaling*res)/np.linalg.norm(scaling*external[columns]));history.append(error)
            if error<1e-9:break
            Kr=K[np.ix_(columns,columns)];dq=np.zeros(len(f));dq[columns]=scaling*np.linalg.solve(Kr*scaling[:,None]*scaling[None,:],-scaling*res)
            merit=U-external@ravel_pose(r);accepted=False
            for bt in range(12):
                dr=dq/2**bt;rr=r.copy();RR=R.copy()
                for i in range(1,n+1):rr[i]+=dr[6*i:6*i+3];RR[i]=R[i]@expm(np.array(b.skew(dr[6*i+3:6*i+6]),float))
                trial=sum(b.evaluate(rr[i:i+2],RR[i:i+2],L,reference_K,np.array([1.,0,0,0,0,0]))['U'] for i in range(n))
                if trial-external@ravel_pose(rr)<=merit+1e-14:r,R=rr,RR;accepted=True;break
            if not accepted:raise ValueError('Newton line search')
        else:raise ValueError('Newton iteration cap')
    return r,R,history
def ravel_pose(r):
    out=np.zeros(6*len(r))
    for i in range(len(r)):out[6*i:6*i+3]=r[i]
    return out
def continuous_static(K,P):
    EA=K[0,0];S=K[1,1];EI=K[5,5]
    def rhs(s,y):
        theta=y[2];v1=1+P*np.sin(theta)/EA;v2=P*np.cos(theta)/S
        dx=np.cos(theta)*v1-np.sin(theta)*v2;dy=np.sin(theta)*v1+np.cos(theta)*v2
        return np.array([dx,dy,y[3]/EI,-P*dx])
    def bc(a,z):return np.array([a[0],a[1],a[2],z[3]])
    x=np.linspace(0,1,41);initial=np.array([x,P*x*x*(3-x)/(6*EI)+P*x/S,P*x*(2-x)/(2*EI),P*(1-x)])
    result=solve_bvp(rhs,bc,x,initial,tol=1e-10,max_nodes=5000)
    if not result.success:raise ValueError('independent BVP')
    return result.y[:,-1],float(np.max(result.rms_residuals))

def main():
    dest=ROOT/'rd/structure/beam/records/beam-se3-physical-v4-result.json';assert not dest.exists();start=time.monotonic();records=[];failures=[]
    try:
        # Analytic constant curvature and twist. Section-energy dual norm is
        # used for force/moment comparison because their units differ.
        cases=[]
        for h in [.1,.01,.0001]:
            K,mass,props=section(h)
            for axis in [1,2]:
                for n,total in [(1,1.8),(4,4.5)]:
                    L=1/n;strain=np.zeros(6);strain[3+axis]=total;d=np.array([1.,0,0,0,0,0])+strain
                    full=L*d;r,Q=b.pose_at(np.zeros(3),np.eye(3),np.array(full,dtype=object),1.);rr=np.array([[0.,0.,0.],b.float_array(r)]);RR=np.array([np.eye(3),b.float_array(Q)]);e=b.evaluate(rr,RR,L,K,np.array([1.,0,0,0,0,0]))
                    expected=K@strain
                    rawforce=np.r_[RR[1].T@e['f'][6:9],e['f'][9:12]];rawdiff=rawforce-expected
                    rawdual=float(np.sqrt(rawdiff@np.linalg.solve(K,rawdiff)/(strain@K@strain)));precision=16
                    if K[0,0]/min(K[4,4],K[5,5])>1e6:
                        angle=Dec(str(total))*Dec(str(L));ss,cc=sincos(angle);k=Dec(str(total));Z=Dec(0);I=Dec(1)
                        if axis==2:
                            pr=[ss/k,(1-cc)/k,Z];qr=[[cc,-ss,Z],[ss,cc,Z],[Z,Z,I]]
                        else:
                            pr=[ss/k,Z,-(1-cc)/k];qr=[[cc,Z,ss],[Z,I,Z],[-ss,Z,cc]]
                        rprec=np.array([[Z,Z,Z],pr],object);qprec=np.array([[[Dec(int(i==j)) for j in range(3)] for i in range(3)],qr],object)
                        pe=hp.evaluate(rprec,qprec,L,K,np.array([1.,0,0,0,0,0]));e={key:np.asarray(pe[key],float) for key in ['f','strain']};precision=80
                    force=np.r_[RR[1].T@e['f'][6:9],e['f'][9:12]];difference=force-expected
                    dual=float(np.sqrt(difference@np.linalg.solve(K,difference)/(strain@K@strain)))
                    exact_error=float(np.max(abs(e['strain']-strain)));assert max(dual,exact_error)<=1e-8
                    # Check a point by an independent sine/cosine circle.
                    rs,Rs=b.pose_at(np.zeros(3),np.eye(3),np.array(d,dtype=object),.63)
                    k=total;s=.63;expected_position=np.array([np.sin(k*s)/k,(1-np.cos(k*s))/k,0.]) if axis==2 else np.array([np.sin(k*s)/k,0.,-(1-np.cos(k*s))/k])
                    geometric=float(np.linalg.norm(b.float_array(rs)-expected_position));assert geometric<1e-10
                    cases.append(dict(h=h,axis=axis,n=n,total_rotation=total,element_rotation=total/n,strain_error=exact_error,wrench_dual_error=dual,circle_error=geometric,precision_digits=precision,binary64_wrench_dual_error=rawdual))
        records.append(dict(case='large_bending_two_axes_three_slenderness',cases=cases));print('PURE_BENDING',max(q['wrench_dual_error'] for q in cases),flush=True)
        # Axial and torsional constant-strain patches.
        for j in [0,3]:
            K,mass,props=section(.1);strain=np.zeros(6);strain[j]=.02;d=np.array([1.,0,0,0,0,0])+strain;r,Q=b.pose_at(np.zeros(3),np.eye(3),d,1.);e=b.evaluate(np.array([np.zeros(3),b.float_array(r)]),np.array([np.eye(3),b.float_array(Q)]),1.,K,np.array([1.,0,0,0,0,0]));assert np.max(abs(e['strain']-strain))<1e-10
        records.append(dict(case='axial_torsion_constant_strain',passed=True))
        linear=[]
        for h in [.1,.01,.0001]:
            K,mass,props=section(h);P=K[5,5]*1e-6;exact=P/(3*K[5,5])+P/K[1,1]
            for n in [4,8,16,32]:
                root,M,_=linear_operators(n,h);F=np.zeros(6*n);F[-5]=P;u,res,cond=linear_solution(root,F);err=abs(u[-5]-exact)/abs(exact)
                linear.append(dict(h=h,n=n,tip=float(u[-5]),exact=exact,error=err,residual=res,root_rcond=cond,state_precision=linear_solution.last_state));assert res<=1e-8
        for h in [.1,.01,.0001]:
            seq=[q['error'] for q in linear if q['h']==h];assert seq[-1]<=.01 and all(a>=c for a,c in zip(seq[:-1],seq[1:]))
        for n in [4,8,16,32]:
            seq=[q['error'] for q in linear if q['n']==n];assert max(seq)/min(seq)<=2
        records.append(dict(case='Timoshenko_tip_force_convergence',cases=linear));print('LINEAR',linear[-1],flush=True)
        # Nonlinear physical tip force, independently integrated continuum ODE.
        K,mass,props=section(.1);P=.36*K[5,5];reference,bvp_residual=continuous_static(K,P);nonlinear_cases=[]
        for n in [4,8,16,32]:
            r,R,history=nonlinear(n,P,K);theta=np.arctan2(R[-1,1,0],R[-1,0,0]);position=np.linalg.norm(r[-1,:2]-reference[:2])/np.linalg.norm(reference[:2]-np.array([1.,0.]));angle=abs(theta-reference[2])/abs(reference[2]);nonlinear_cases.append(dict(n=n,position_error=float(position),angle_error=float(angle),residual=history[-1],iterations=len(history)))
            assert history[-1]<=1e-8;print('NONLINEAR',nonlinear_cases[-1],flush=True)
        assert max(nonlinear_cases[-1][key] for key in ['position_error','angle_error'])<=.01
        for key in ['position_error','angle_error']:assert all(a[key]>=c[key] for a,c in zip(nonlinear_cases[:-1],nonlinear_cases[1:]))
        records.append(dict(case='nonlinear_tip_force_vs_continuum_BVP',cases=nonlinear_cases,continuum_tip=reference.tolist(),continuum_residual=bvp_residual))
        # Independent analytic fixed-free rods and continuum Timoshenko modes.
        modes=[];h=.001;K,mass,props=section(h);bendref=continuous_bending(props,h)
        for n in [16,32,64]:
            root,M,_=linear_operators(n,h,True)
            for name,indices,c in [('axial',[6*i for i in range(n)],np.sqrt(K[0,0]/mass[0,0])),('torsion',[6*i+3 for i in range(n)],np.sqrt(K[3,3]/mass[3,3]))]:
                omega,back,orth=frequencies(root,M,np.array(indices));exact=np.arange(1,6,2)*np.pi/2*c;errors=abs(omega[:3]-exact)/exact
                modes.append(dict(n=n,type=name,frequency_errors=errors.tolist(),backward=back,mass_orthogonality=orth));assert max(back,orth)<=1e-8
                if n>=32:assert max(errors)<=.004
            indices=np.array([6*i+j for i in range(n) for j in [1,5]]);omega,back,orth=frequencies(root,M,indices);errors=abs(omega[:2]-bendref)/bendref
            modes.append(dict(n=n,type='bending',frequency_errors=errors.tolist(),backward=back,mass_orthogonality=orth));assert max(back,orth)<=1e-8
            if n==64:assert max(errors)<=.004
        records.append(dict(case='consistent_modes_vs_analytic_and_continuum',cases=modes,bending_reference=bendref.tolist()));print('MODES',modes[-1],flush=True)
        # Independent Saint-Venant cross-section Poisson problem.
        width=.2;height=.1;exact=torsion_rectangle(width,height,401);torsion=[]
        for N in [32,64,128]:
            nx=2*N;ny=N;dx=width/nx;dy=height/ny;Tx=diags([-np.ones(nx-2),2*np.ones(nx-1),-np.ones(nx-2)],[-1,0,1])/dx**2;Ty=diags([-np.ones(ny-2),2*np.ones(ny-1),-np.ones(ny-2)],[-1,0,1])/dy**2
            A=kron(eye(ny-1),Tx)+kron(Ty,eye(nx-1));psi=spsolve(A.tocsc(),np.full((nx-1)*(ny-1),2.));J=2*sum(psi)*dx*dy;torsion.append(dict(N=N,J=float(J),relative_error=float(abs(J-exact)/exact)))
        assert torsion[-1]['relative_error']<=.01 and all(a['relative_error']>=c['relative_error'] for a,c in zip(torsion[:-1],torsion[1:]));polar=width*height*(width**2+height**2)/12;negative=abs(polar-exact)/exact;assert negative>.1
        records.append(dict(case='Saint_Venant_torsion_vs_2D_Poisson',cases=torsion,series_J=exact,polar_substitution_defect=negative))
    except Exception as exc:failures.append(dict(reason=repr(exc),traceback=traceback.format_exc()));print('FAIL',failures[-1],flush=True)
    paths=[Path(__file__),ROOT/'rd/structure/beam/experiments/check_beam_se3_physical_v3.py',ROOT/'rd/structure/beam/records/beam-se3-physical-v3-result.json',ROOT/'rd/structure/beam/experiments/check_beam_se3_physical_v2.py',ROOT/'rd/structure/beam/records/beam-se3-physical-v2-result.json',Path(hp.__file__),Path(hp.ad.__file__),ROOT/'rd/structure/beam/experiments/beam_decimal_reference.py',ROOT/'rd/structure/beam/records/beam-se3-precision-addendum.md',ROOT/'rd/structure/beam/experiments/check_beam_se3_physical.py',ROOT/'rd/structure/beam/records/beam-se3-physical-result.json',Path(b.__file__),Path(b.ad.__file__),ROOT/'rd/structure/beam/records/beam-se3-adoption-protocol.md']
    result=dict(passed=len(records)==6 and not failures,records=records,failures=failures,seconds=time.monotonic()-start,
        sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print('BEAM_PHYSICAL_RESULT',result['passed'],flush=True)
if __name__=='__main__':main()
