import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import time,json,hashlib,numpy as np
from scipy.linalg import solve,lu_factor,lu_solve,eig
from numba import njit,prange
from numpy.polynomial.legendre import leggauss
from uvlm_frequency_fast import geometry,bound,age_nodes
from uvlm_age_scalar import point_derivative,age_integral
ROOT=rd_paths.ROOT
@njit(parallel=True)
def nodal_wake(P,N,te,eps,h,R,gx,gw):
    n=len(te)-1;W=np.zeros((len(P),R+1,n))
    for k in prange(len(P)):
        for j in range(n):
            dl=te[j+1]-te[j];L=np.sqrt(dl@dl);e=dl/L;r=P[k]-te[j]
            for row in range(R):
                for z,w in zip(gx,gw):
                    s=(z+1)/2;age=(row+s)*h;v=20*h*w/2*point_derivative(r[0]-20*age,r[1],r[2],e[0],e[1],e[2],L,eps,N[k,0],N[k,1],N[k,2]);W[k,row,j]+=(1-s)*v;W[k,row+1,j]+=s*v
    return W
def operators():
    g=geometry(2,4);X=g['X'];A=bound(g['P'],g['normal'],X,1e-5);nb=len(A);C=np.zeros((2,nb));H=np.zeros((2,nb));n=g['n']
    for j in range(nb):
        gamma=np.zeros(nb);gamma[j]=1;G=gamma.reshape(g['m'],n);horse=(G-np.vstack([np.zeros((1,n)),G[:-1]])).ravel();f=1.2*horse[:,None]*np.cross([20.,0,0],g['right']-g['left']);C[:,j]=[f[:,2].sum(),np.cross((g['left']+g['right'])/2-[.375,0,0],f)[:,1].sum()]
    quads=np.stack([X[:-1,:-1],X[:-1,1:],X[1:,1:],X[1:,:-1]],axis=2).reshape(-1,4,3)
    for ids in [[0,2,1],[0,3,2]]:
        tri=quads[:,ids];f=1.2*np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])/2;H[0]+=f[:,2];H[1]+=np.cross(tri.mean(axis=1)-[.375,0,0],f)[:,1]
    Bq=np.column_stack([np.zeros(nb),-20*g['normal'][:,2]]);Bv=g['normal'][:,2,None]*np.column_stack([np.ones(nb),-(g['P'][:,0]-.375)]);return g,A,C,H,Bq,Bv
def main():
    start=time.monotonic();g,A,C,H,Bq,Bv=operators();nb=len(A);n=g['n'];M=np.diag([40.,8.]);K=np.diag([6000.,20000.]);D=np.diag([10.,20.]);omega=10.;f=np.array([1.,.2]);refs=[]
    for order in [16,24]:
        ages,weights=age_nodes(g,omega,20.,order);W=age_integral(g['P'],g['normal'],g['X'][-1],1e-5,ages,weights);AC=A.astype(complex);AC[:,-n:]+=W;L=solve(AC,Bq+1j*omega*Bv);q=solve(-omega**2*M+1j*omega*D+K-(C+1j*omega*H)@L,f);gamma=L@q;refs.append(np.r_[q,1j*omega*q,(C+1j*omega*H)@gamma])
    quadrature=float(np.max(abs((refs[0]-refs[1])/refs[1])));assert quadrature<1e-9;records=[];errors=[];gx,gw=leggauss(12)
    for R in [32,64,128,256]:
        h=1/R;alpha=1.5/h;W=nodal_wake(g['P'],g['normal'],g['X'][-1],1e-5,h,R,gx,gw);z=np.exp(1j*omega*h);delta=(3-4/z+1/z**2)/(2*h);Wh=np.einsum('prn,r->pn',W,z**(-np.arange(R+1)));AC=A.astype(complex);AC[:,-n:]+=Wh;L=solve(AC,Bq+delta*Bv);q=solve(delta**2*M+delta*D+K-(C+delta*H)@L,f);gamma=L@q;obs=np.r_[q,delta*q,(C+delta*H)@gamma];error=abs((obs-refs[-1])/refs[-1]);errors.append(error)
        S=np.block([[alpha*np.eye(2),-np.eye(2),np.zeros((2,nb))],[K,alpha*M+D,-C-alpha*H],[-Bq,-Bv,A+np.pad(W[:,0],((0,0),(nb-n,0)))]]);factor=lu_factor(S);size=4+nb
        # State [q,v,Gamma] now and previous, plus free wake nodes 1..R.
        def step(state,forcing):
            old=state[:size];prev=state[size:2*size];mu=state[2*size:].reshape(R,n,-1);newmu=np.concatenate([old[-n:][None],mu[:-1]],axis=0);history=(-4*old+prev)/(2*h);rhs=np.concatenate([-history[:2],forcing-M@history[2:4]+H@history[4:],-np.einsum('prn,rnk->pk',W[:,1:],newmu)],axis=0);new=lu_solve(factor,rhs);return np.concatenate([new,old,newmu.reshape(R*n,-1)],axis=0)
        one=np.r_[q,delta*q,gamma];state=np.concatenate([one,one/z,(z**(-np.arange(1,R+1)))[:,None]*gamma[-n:][None]],axis=None)[:,None];initial=state.copy()
        for j in range(1,21):state=step(state,(f*z**j)[:,None])
        time_error=float(np.linalg.norm(state[:size,0]-one*z**20)/np.linalg.norm(one));assert time_error<1e-9
        resumed=initial.copy()
        for j in range(1,11):resumed=step(resumed,(f*z**j)[:,None])
        restored=resumed.copy()
        for j in range(11,21):restored=step(restored,(f*z**j)[:,None])
        assert np.array_equal(state,restored);row=dict(intervals=R,relative_errors=error.tolist(),time_frequency_difference=time_error,restart_bitwise=True)
        if R<=64:
            dim=len(state);Smap=step(np.eye(dim),np.zeros((2,dim)));rng=np.random.default_rng(R);u=rng.normal(size=(dim,1));v=rng.normal(size=(dim,1));hh=1e-5;fd=(step(u+hh*v,np.zeros((2,1)))-step(u-hh*v,np.zeros((2,1))))/(2*hh);jerr=float(np.linalg.norm(fd-Smap@v)/np.linalg.norm(Smap@v));assert jerr<1e-8;ev,V=eig(Smap);pick=int(np.argmax(abs(ev)));mode=V[:,pick:pick+1];advanced=mode.copy()
            for _ in range(12):advanced=step(advanced,np.zeros((2,1)))
            eigerr=float(np.linalg.norm(advanced-ev[pick]**12*mode)/np.linalg.norm(advanced));assert eigerr<1e-9;row.update(full_step_dimension=dim,step_JVP_difference=jerr,eigenmode_time_difference=eigerr,spectral_radius=float(max(abs(ev))))
        records.append(row);assert time.monotonic()-start<180
    orders=np.log2(np.array(errors[-2])/np.array(errors[-1]));passed=bool(np.all((orders>1.8)&(orders<2.2)) and max(errors[-1])<.01);dest=ROOT/'rd/coupled/j1/records/coupled-bdf2-prescribed-result.json';assert not dest.exists();paths=[Path(__file__),ROOT/'rd/coupled/j1/records/coupled-bdf2-prescribed-protocol.md',ROOT/'rd/coupled/j1/experiments/uvlm_age_scalar.py',ROOT/'rd/aero/unsteady-reference/experiments/uvlm_frequency_fast.py'];result=dict(passed=passed,scope='linear heave/pitch with prescribed material wake, fixed physical tail; not general J1 adoption',quadrature_difference=quadrature,records=records,final_orders=orders.tolist(),seconds=time.monotonic()-start,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths});dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2));assert passed
if __name__=='__main__':main()
