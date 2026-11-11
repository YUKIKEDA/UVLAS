import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import json,hashlib,numpy as np
from scipy.spatial.transform import Rotation
ROOT=rd_paths.ROOT
def skew(v):
    x,y,z=v;return np.array([[0,-z,y],[z,0,-x],[-y,x,0]],dtype=np.result_type(v,float))
def left(v):
    a=np.sqrt(np.sum(v*v));S=skew(v)
    if abs(a)<1e-5:A=.5-a*a/24+a**4/720;B=1/6-a*a/120+a**4/5040
    else:A=(1-np.cos(a))/(a*a);B=(a-np.sin(a))/(a*a*a)
    return np.eye(3)+A*S+B*S@S
def left_dot(v,d):return left(v.astype(complex)+1j*1e-25*d).imag/1e-25
def path(t):
    v=np.array([.3*np.sin(t),.2*np.cos(.7*t),.1*np.sin(.3*t)]);d=np.array([.3*np.cos(t),-.14*np.sin(.7*t),.03*np.cos(.3*t)]);dd=np.array([-.3*np.sin(t),-.098*np.cos(.7*t),-.009*np.sin(.3*t)]);return v,d,dd
def chart(R,R0,w):
    if np.linalg.norm(R.T@R-np.eye(3))>1e-12 or np.linalg.det(R)<0:raise ValueError('invalid_rotation')
    v=Rotation.from_matrix(R0.T@R).as_rotvec()
    if np.linalg.norm(v)>np.pi-.02:raise ValueError('log_branch_margin')
    return v,np.linalg.solve(left(v),R0.T@w)
def main():
    R0=Rotation.from_rotvec([.1,.3,-.2]).as_matrix();R1=R0@Rotation.from_rotvec([.6,-.2,.1]).as_matrix();records=[];t=.8;theta,d,dd=path(t);R=R0@Rotation.from_rotvec(theta).as_matrix();w=R0@left(theta)@d;wdot=R0@(left_dot(theta,d)@d+left(theta)@dd);tt,vv=chart(R,R1,w);pose_error=np.linalg.norm(R1@Rotation.from_rotvec(tt).as_matrix()-R);velocity_error=np.linalg.norm(R1@left(tt)@vv-w);assert max(pose_error,velocity_error)<1e-12
    inertia=R@np.diag([2.,3.,4.])@R.T;old=.5*d@(left(theta).T@R0.T@inertia@R0@left(theta))@d;new=.5*vv@(left(tt).T@R1.T@inertia@R1@left(tt))@vv;assert abs(old-new)<1e-12;moment=np.array([.2,-.1,.3]);power_old=d@(left(theta).T@R0.T@moment);power_new=vv@(left(tt).T@R1.T@moment);assert abs(power_old-power_new)<1e-12;records.append(dict(case='physical_pose_twist_energy_and_power',pose_error=float(pose_error),twist_error=float(velocity_error),energy_difference=float(abs(old-new)),power_difference=float(abs(power_old-power_new))))
    errors=[];wrong=None
    for h in [.1,.05,.025,.0125]:
        olds=[];news=[];angles=[]
        for tau in [t,t-h,t-2*h]:
            a,b,_=path(tau);rr=R0@Rotation.from_rotvec(a).as_matrix();ww=R0@left(a)@b;aa,bb=chart(rr,R1,ww);olds.append(b);news.append(bb);angles.append(aa)
        ao=(3*olds[0]-4*olds[1]+olds[2])/(2*h);an=(3*news[0]-4*news[1]+news[2])/(2*h);physical_old=R0@(left_dot(theta,olds[0])@olds[0]+left(theta)@ao);physical_new=R1@(left_dot(angles[0],news[0])@news[0]+left(angles[0])@an);error=max(np.linalg.norm(physical_old-wdot),np.linalg.norm(physical_new-wdot));errors.append(error)
        aw=(3*news[0]-4*olds[1]+olds[2])/(2*h);wrong=R1@(left_dot(angles[0],news[0])@news[0]+left(angles[0])@aw);records.append(dict(case='both_histories_remapped',h=h,absolute_acceleration_error=float(error),old_new_algorithmic_difference=float(np.linalg.norm(physical_old-physical_new))))
    order=np.log2(errors[-2]/errors[-1]);assert 1.8<order<2.2 and errors[-1]/np.linalg.norm(wdot)<.01;negative=float(np.linalg.norm(wrong-wdot)/np.linalg.norm(wdot));assert negative>.1;rejected=[]
    for label,rr in [('branch',R0@Rotation.from_rotvec([np.pi-.01,0,0]).as_matrix()),('improper',np.diag([-1.,1.,1.]))]:
        try:chart(rr,R0,np.zeros(3))
        except ValueError:rejected.append(label)
        else:raise AssertionError(label)
    records.append(dict(case='current_only_recenter_negative_control',wrong_relative_error=negative,rejections=rejected));paths=[Path(__file__),ROOT/'docs/specs/coupled/j1-discrete-equilibrium-and-chart-contract.md'];result=dict(passed=True,records=records,final_order=float(order),scope='SO3 chart component of SE3 pose/twist histories, physical inertia/power and branch guards',sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths});dest=ROOT/'rd/coupled/j1/records/so3-time-history-result.json';assert not dest.exists();dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
