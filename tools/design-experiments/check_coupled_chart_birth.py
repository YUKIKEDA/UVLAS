from pathlib import Path
import json,hashlib,time,copy,numpy as np
from scipy.integrate import solve_ivp
from run_frequency_delta import peak_bytes
ROOT=Path(__file__).resolve().parents[2]
def pivot(t):return np.array([.1*np.sin(t),.07*np.cos(.8*t)])
def pdot(t):return np.array([.1*np.cos(t),-.056*np.sin(.8*t)])
def pddot(t):return np.array([-.1*np.sin(t),-.0448*np.cos(.8*t)])
def radial(r):return np.array([np.cos(r),np.sin(r)])
def tangent(r):return np.array([-np.sin(r),np.cos(r)])
def torque(r,u,t):return -np.sin(r)-.2*u+np.cos(.7*t)
def config(r,t):return pivot(t)+radial(r)
def velocity(r,u,t):return pdot(t)+tangent(r)*u
def fluid(x,r,t):return np.array([1+.2*np.sin(x[1])+.1*np.cos(r),.3*np.sin(x[0])+.1*np.sin(r)])
def fluid_jac(x):return np.array([[0,.2*np.cos(x[1])],[.3*np.cos(x[0]),0]])
def reference(tol):
    sol=solve_ivp(lambda t,y:[y[1],torque(y[0],y[1],t)],[0,1],[.3,.2],method='DOP853',rtol=tol,atol=tol*.01,dense_output=True);points=[]
    for birth in [0,.25,.5,.75]:
        initial=config(sol.sol(birth)[0],birth);x=solve_ivp(lambda t,x:fluid(x,sol.sol(t)[0],t),[birth,1],initial,method='DOP853',rtol=tol,atol=tol*.01);points.append(x.y[:,-1])
    r,u=sol.y[:,-1];return np.r_[config(r,1),velocity(r,u,1),np.array(points).ravel()]
def initial(h):
    r,u=.3,.2;v=velocity(r,u,0);a=pddot(0)+tangent(r)*torque(r,u,0)-radial(r)*u*u;return dict(r=r,u=u,rprev=r-h*u,v=v,vprev=v-h*a,points={},index=0)
def step(parent,h,N):
    s=copy.deepcopy(parent);i=parent['index'];t0=i*h;t1=(i+1)*h;alpha=1.5/h;hist=(-4*parent['r']+parent['rprev'])/(2*h);vhist=(-4*parent['v']+parent['vprev'])/(2*h)
    def residual(r):
        u=alpha*r+hist;v=velocity(r,u,t1);return tangent(r)@(alpha*v+vhist-pddot(t1))-torque(r,u,t1)
    r=parent['r']+h*parent['u']
    for k in range(12):
        R=residual(r);der=residual(r+1j*1e-25).imag/1e-25;increment=-R/der;r+=increment
        if abs(increment)<1e-14:break
    else:raise ValueError('mechanical_Newton_failed')
    u=alpha*r+hist;v=velocity(r,u,t1);s.update(r=float(r),u=float(u),rprev=parent['r'],v=v,vprev=parent['v'],index=i+1)
    if i in [0,N//4,N//2,3*N//4]:
        key=str(i*h)
        if key in s['points']:raise ValueError('duplicate_birth')
        s['points'][key]=dict(x=config(parent['r'],t0),previous=None,newborn=True)
    for key,p in s['points'].items():
        old=p['x'].copy();guess=old+h*fluid(old,parent['r'],t0)
        if p['newborn']:
            oldV=fluid(old,parent['r'],t0)
            def f(x):return x-old-h/2*(oldV+fluid(x,r,t1))
            def J(x):return np.eye(2)-h/2*fluid_jac(x)
        else:
            if p['previous'] is None:raise ValueError('missing_material_history')
            history=(-4*old+p['previous'])/(2*h)
            def f(x):return alpha*x+history-fluid(x,r,t1)
            def J(x):return alpha*np.eye(2)-fluid_jac(x)
        for _ in range(12):
            inc=np.linalg.solve(J(guess),-f(guess));guess+=inc
            if np.linalg.norm(inc)<1e-14:break
        else:raise ValueError('material_Newton_failed')
        p.update(x=guess,previous=old,newborn=False)
    q=config(r,t1);G=radial(r);cres=abs(np.linalg.norm(q-pivot(t1))-1);vres=abs(G@(v-pdot(t1)));assert max(cres,vres)<1e-12;return s
def output(s):return np.r_[config(s['r'],1),s['v'],np.array([p['x'] for p in s['points'].values()]).ravel()]
def main():
    start=time.monotonic();ref=reference(5e-13);coarse_ref=reference(2e-12);refdiff=np.linalg.norm(ref-coarse_ref);assert refdiff<1e-9;records=[];errors=[]
    for N in [20,40,80,160]:
        h=1/N;s=initial(h)
        for i in range(N):
            s=step(s,h,N)
            if i==N//2-1:saved=copy.deepcopy(s)
        resumed=copy.deepcopy(saved)
        for i in range(N//2,N):resumed=step(resumed,h,N)
        assert np.array_equal(output(s),output(resumed));value=output(s);groups=[np.linalg.norm(value[:2]-ref[:2]),np.linalg.norm(value[2:4]-ref[2:4]),np.max(np.linalg.norm((value[4:]-ref[4:]).reshape(4,2),axis=1))];errors.append(groups);bad=copy.deepcopy(saved);next(iter(bad['points'].values()))['previous']=None
        try:step(bad,h,N)
        except ValueError:pass
        else:raise AssertionError('missing_material_history')
        records.append(dict(steps=N,absolute_errors=groups,restart_bitwise=True,missing_material_history_rejected=True))
    orders=np.log2(np.array(errors[-2])/np.array(errors[-1]));peak=peak_bytes();passed=bool(min(orders)>1.8 and max(orders)<2.2 and max(errors[-1])<.01 and peak<2*2**30 and time.monotonic()-start<180);paths=[Path(__file__),ROOT/'docs/verification/coupled-chart-birth-protocol.md',ROOT/'docs/coupled-time-stability-specification.md'];result=dict(passed=passed,scope='manufactured nonlinear moving constraint chart and solenoidal birth/transport; not nonlinear UVLM physical validation',records=records,final_orders=orders.tolist(),reference_sensitivity=float(refdiff),seconds=time.monotonic()-start,peak_bytes=peak,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths});dest=ROOT/'docs/verification/coupled-chart-birth-result.json';assert not dest.exists();dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2));assert passed
if __name__=='__main__':main()
