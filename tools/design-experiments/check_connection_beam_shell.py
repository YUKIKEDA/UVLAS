"""Actual nonlinear beam-shell KKT versus independent nonlinear elimination."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,time,traceback
import numpy as np
from scipy.linalg import solve,null_space,block_diag,eigh,qr,cholesky,solve_triangular,svd
from scipy.optimize import brentq
import connection_shell as pg
import connection_exact_trace as exact
import connection_kinematics as kin
import beam_se3_v2 as beam
ROOT=Path(__file__).resolve().parents[2]

def section():
    A=.05;I=A*A/(4*np.pi);K=np.diag([1000*A,500*A*5/6,500*A*5/6,1000*I,1000*I,1000*I]);M=beam.spatial_mass(2*A,np.zeros(3),np.diag([4*I,2*I,2*I]));return K,M

def prepare():
    m=pg.flat_model(1,.1,1.,1.);d=pg.data(m,origin=np.array([1.,-.5,0.]));ys=[F(2*j+1,40) for j in range(20)];points=[(F(0),y,z) for y in ys for z in [F(-1),F(0),F(1)]];T,record=exact.assemble([m],[[(0,points,1)]]);Z=np.kron(null_space(T),np.eye(3));p0=np.array([1.,0.,0.]);eta=d['q0']-d['translation'][:,None]*p0
    p,w,area=pg.face(m,axis=0,side=1.);N,_,_=pg.point(m,p);face=N.T@(w*np.linalg.norm(area,axis=1)/.1)
    B=np.einsum('kab,pib->pkia',pg.material.MB,d['G'],optimize=True).reshape(len(d['G']),6,-1);root=(np.sqrt(d['weights']*1000)[:,None,None]*B).reshape(-1,3*m['ns']);Rroot=qr(root,mode='economic',check_finite=False)[1]
    return dict(d=d,T=T,Z=Z,p0=p0,eta=eta,weights=T@d['translation'],face=face,trace=record,root=Rroot)

def constraints(data,p,R,u):
    T=data['T'];eta=T@data['eta'];weights=data['weights'];nr=len(T);C=T@u-weights[:,None]*(p-data['p0'])-(eta@R.T-eta);J=np.zeros((3*nr,6+u.size));J[:,6:]=np.kron(T,np.eye(3));H=np.zeros((3*nr,3,3))
    for i in range(nr):
        J[3*i:3*i+3,:3]=-weights[i]*np.eye(3);J[3*i:3*i+3,3:6]=R@np.array(beam.skew(eta[i]),float)
        for a in range(3):
            A=np.array(beam.skew(np.eye(3)[a]),float)
            for b in range(3):
                B=np.array(beam.skew(np.eye(3)[b]),float);H[3*i:3*i+3,a,b]=-R@((A@B+B@A)/2)@eta[i]
    return C.ravel(),J,H

def internal(data,p,R,u):
    kb,_=section();e=beam.evaluate(np.array([[0.,0.,0.],p]),np.array([np.eye(3),R]),1.,kb,np.array([1.,0,0,0,0,0]));U,fs,Ks,info=pg.material.evaluate(data['d'],u,True);assert info['min_det']>0
    return e['U']+U,np.r_[e['f'][6:],fs.ravel()],block_diag(e['K'][6:,6:],Ks)

def update(p,R,dq):return p+dq[:3],R@beam.float_array(beam.exp_rotation(dq[3:6]))

def solve_case(data,F,method):
    Z=data['Z'];d=data['d'];p=data['p0']+np.array([.0001,-.00002,.00003]);R=beam.float_array(beam.exp_rotation(np.array([.003,-.002,.001])));v=np.zeros(Z.shape[1]);value,_,_=kin.shoe(d['q0'],d['translation'],data['p0'],p,R);u=value-d['q0'];lam=np.zeros(3*len(data['T']));external=np.r_[np.zeros(6),(data['face'][:,None]*F).ravel()];history=[]
    for fraction in [.25,.5,1.]:
        for iteration in range(20):
            if method=='elimination':
                value,Gs,Hs=kin.shoe(d['q0'],d['translation'],data['p0'],p,R);u=value-d['q0']+(Z@v).reshape(-1,3);G=np.zeros((6+u.size,6+len(v)));G[:6,:6]=np.eye(6);G[6:,:6]=Gs;G[6:,6:]=Z
            U,f,K=internal(data,p,R,u);res=f-fraction*external;C,J,H=constraints(data,p,R,u)
            if method=='kkt':
                lag=res+J.T@lam;KL=K.copy();KL[3:6,3:6]+=np.einsum('i,iab->ab',lam,H);rr=np.r_[lag,C];A=np.block([[KL,J.T],[J,np.zeros((len(C),len(C)))]])
                err=float(np.linalg.norm(lag)/max(np.linalg.norm(fraction*external),1e-30));gap=float(np.max(abs(C)))
            else:
                rr=G.T@res;A=G.T@K@G;A[:6,:6]+=np.einsum('i,iab->ab',res[6:],Hs);err=float(np.linalg.norm(rr)/max(np.linalg.norm(fraction*external),1e-30));gap=float(np.max(abs(C)))
            history.append(dict(fraction=fraction,iteration=iteration,residual=err,gap=gap));print(method,fraction,iteration,err,gap,flush=True)
            if err<2e-10 and gap<1e-11:break
            scales=1/np.sqrt(np.maximum(abs(np.diag(A)),1e-12));dx=scales*solve(A*scales[:,None]*scales[None,:],-rr*scales,assume_a='sym');p,R=update(p,R,dx[:6])
            if method=='kkt':u=u+dx[6:6+u.size].reshape(-1,3);lam+=dx[6+u.size:]
            else:v+=dx[6:]
        else:raise ValueError('nonlinear_equilibrium_not_converged')
    if method=='elimination':lam=np.linalg.lstsq(J.T,-res,rcond=None)[0]
    raw=float(np.linalg.norm(res+J.T@lam)/np.linalg.norm(external));assert raw<1e-8 and np.max(abs(C))<1e-9
    return dict(p=p,R=R,u=u,lam=lam,raw=raw,history=history,external=external,K=K,J=J)

def modal(data,n):
    kb,mb=section();L=1/n;r=np.array([[0.,0.,0.],[L,0,0]]);Rs=np.array([np.eye(3)]*2);local=beam.evaluate(r,Rs,L,kb,np.array([1.,0,0,0,0,0]))['K'];me,_=beam.inertia(r,Rs,L,mb,6);K=np.zeros((6*(n+1),6*(n+1)));M=K.copy()
    for i in range(n):K[6*i:6*i+12,6*i:6*i+12]+=local;M[6*i:6*i+12,6*i:6*i+12]+=me
    localB=beam.evaluate(r,Rs,L,kb,np.array([1.,0,0,0,0,0]))['B'];rb=np.zeros((6*n,6*(n+1)))
    for i in range(n):rb[6*i:6*i+6,6*i:6*i+12]=np.sqrt(L)*np.linalg.cholesky(kb).T@localB
    root=block_diag(rb[:,6:],data['root']);K=block_diag(K[6:,6:],pg.linear_stiffness(data['d']));M=block_diag(M[6:,6:],np.kron(data['d']['M'],np.eye(3)));_,J0,_=constraints(data,data['p0'],np.eye(3),np.zeros_like(data['d']['q0']));J=np.zeros((len(J0),len(K)));J[:,-(len(J0[0])):]=J0;Z=null_space(J);G=np.zeros((len(K),6*n+data['Z'].shape[1]));G[:6*n,:6*n]=np.eye(6*n);_,Gs,_=kin.shoe(data['d']['q0'],data['d']['translation'],data['p0'],data['p0'],np.eye(3));G[6*n:,6*n-6:6*n]=Gs;G[6*n:,6*n:]=data['Z']
    def factored_modes(A):
        Lm=cholesky(A.T@M@A,lower=True);B=solve_triangular(Lm,(root@A).T,lower=True).T;_,s,Vt=svd(B,full_matrices=False);V=solve_triangular(Lm.T,Vt.T[:,::-1],lower=False);return s[::-1]**2,A@V
    w,modes=factored_modes(Z);w2,_=factored_modes(G);spectral=float(np.max(abs(w-w2)/np.maximum(abs(w),1.)));orth=float(np.linalg.norm(modes.T@M@modes-np.eye(len(w)),np.inf));residual=float(np.linalg.norm(Z.T@(K@modes-M@modes*w))/np.linalg.norm(K@modes));angle=np.arctan(1/np.sqrt(2));ks=[angle,np.pi-angle,np.pi+angle];errors=[];macs=[];indices=[]
    for k in ks:
        target=np.zeros(len(K));x=np.arange(1,n+1)/n;target[:6*n:6]=np.sin(k*x)/(50*k);xx=data['d']['p'][:,0];a=np.sin(k)/(50*k);N0=np.cos(k);values=a*np.cos(k*xx)+N0*np.sin(k*xx)/(100*k);q=np.linalg.solve(data['d']['M'],data['d']['N'].T@(2*data['d']['weights']*values));target[6*n::3]=q;mac=(modes.T@M@target)**2/(target@M@target);j=int(np.argmax(mac));indices.append(j);macs.append(float(mac[j]));errors.append(float(abs(np.sqrt(w[j])/(np.sqrt(500)*k)-1)))
    print('MODAL',n,dict(errors=errors,macs=macs,indices=indices,spectral=spectral,orth=orth,residual=residual),flush=True)
    assert spectral<1e-8 and orth<1e-8 and residual<1e-8 and len(set(indices))==3 and min(macs)>.9
    return dict(n=n,axial_frequency_errors=errors,mac=macs,mode_indices=indices,spectrum_difference=spectral,mass_orthogonality=orth,original_reduced_eigen_residual=residual)

def main():
    start=time.monotonic();records=[];failures=[];data=None
    try:
        data=prepare()
        for F in [np.array([1.,0.,0.]),np.array([1.,.004,-.003])]:
            a=solve_case(data,F,'kkt');b=solve_case(data,F,'elimination');difference=float(max(np.max(abs(a[k]-b[k])) for k in ['p','R','u']));dual=float(np.linalg.norm(a['J'].T@(a['lam']-b['lam']))/np.linalg.norm(a['external']));assert difference<1e-8 and dual<1e-8
            row=dict(force=F.tolist(),kkt_residual=a['raw'],elimination_residual=b['raw'],physical_state_difference=difference,reaction_difference=dual,kkt_history=a['history'],elimination_history=b['history'])
            if F[1]==0:
                lam=brentq(lambda x:50*x*(x*x-1)-1,1,1.1);tipx=data['d']['N']@a['u'];p,w,area=pg.face(data['d']['model'],axis=0,side=1.);N,_,_=pg.point(data['d']['model'],p);tip=float(np.mean((N@a['u'])[:,0]));expected=1/50+lam-1;row['series_tip']=tip;row['independent_series_tip']=expected;row['series_error']=abs(tip/expected-1);assert row['series_error']<1e-8
            records.append(row)
        modalrows=[modal(data,n) for n in [4,8,16]];assert max(modalrows[-1]['axial_frequency_errors'])<.004 and min(modalrows[-1]['mac'])>.99;records.append(dict(modal=modalrows))
    except Exception as exc:failures.append(dict(reason=repr(exc),traceback=traceback.format_exc()));print('FAIL',failures[-1],flush=True)
    # Publication is separate from exploratory execution; failed development runs
    # are retained as logs, never presented as acceptance evidence.
    import sys
    if '--publish' in sys.argv:
        paths=[Path(__file__),Path(pg.__file__),Path(exact.__file__),Path(kin.__file__),Path(beam.__file__),ROOT/'docs/verification/connection-adoption-protocol.md'];dest=ROOT/'docs/verification/connection-beam-shell-result.json';assert not dest.exists();result=dict(passed=not failures and len(records)==3,records=records,failures=failures,seconds=time.monotonic()-start,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths});dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print('CONNECTION_BEAM_SHELL_RESULT',result['passed'],flush=True)
    else:print(json.dumps(dict(records=records,failures=failures),indent=2),flush=True)
if __name__=='__main__':main()
