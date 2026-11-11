"""New element-independent Mandel, reference axis, Schur and stack identities."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json
import numpy as np
ROOT=rd_paths.ROOT
MB=np.zeros((6,3,3))
for i in range(3):MB[i,i,i]=1.
for k,(i,j) in enumerate([(1,2),(0,2),(0,1)],3):MB[k,i,j]=MB[k,j,i]=1/np.sqrt(2.)
D=np.diag([1.,1.,1.,np.sqrt(2.),np.sqrt(2.),np.sqrt(2.)])
def rel(a,b):return float(np.linalg.norm(a)/max(1.,np.linalg.norm(b)))
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def tensor(e):return np.einsum('k,kij->ij',e,MB)
def mandel(E):return np.einsum('kij,ij->k',MB,E)
def axis(J,w):
    cross=np.cross(J[:,0],J[:,1]);norm=np.linalg.norm(cross)
    if norm==0 or np.linalg.det(J)<=0:raise ValueError('invalid_reference_frame')
    n=cross/norm;v=w-n*(w@n)
    if np.linalg.norm(v)==0:raise ValueError('invalid_fiber_projection')
    a=v/np.linalg.norm(v)
    return np.column_stack((a,np.cross(n,a),n))
def transform(A):return np.einsum('kij,ia,lab,jb->kl',MB,A,MB,A)
def rotation(angle):return np.array([[np.cos(angle),-np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1.]])
def moments(layers):
    A=np.zeros((6,6));B=A.copy();DD=A.copy();N=np.zeros(6);M=N.copy();constant=0.
    for lo,hi,Ceng,R,e0 in layers:
        T=transform(R);C=T@D@Ceng@D@T.T;e=T@e0;c=C@e
        t=hi-lo;first=(hi*hi-lo*lo)/2;second=(hi**3-lo**3)/3
        A+=t*C;B+=first*C;DD+=second*C;N+=t*c;M+=first*c;constant+=t*(e@c)/2
    return np.block([[A,B],[B,DD]]),np.r_[N,M],constant
def direct_stack(layers,a,b):
    energy=0.;force=np.zeros(12);points,weights=np.polynomial.legendre.leggauss(2)
    for lo,hi,Ceng,R,e0 in layers:
        for p,w in zip(points,weights):
            z=(lo+hi)/2+(hi-lo)*p/2;weight=w*(hi-lo)/2
            E=R.T@tensor(a+z*b)@R-tensor(e0)
            eng=np.array([E[0,0],E[1,1],E[2,2],2*E[1,2],2*E[0,2],2*E[0,1]])
            stress=Ceng@eng;S=np.array([[stress[0],stress[5],stress[4]],[stress[5],stress[1],stress[3]],[stress[4],stress[3],stress[2]]])
            s=mandel(R@S@R.T);energy+=weight*(eng@stress)/2;force+=weight*np.r_[s,z*s]
    return energy,force
def main():
    dest=ROOT/'rd/structure/material/records/material-layer-kernel-result.json';assert not dest.exists();records=[]
    H=np.diag([1/130,1/10,1/8,1/3.5,1/4.5,1/5]);H[0,1]=H[1,0]=-.28/130;H[0,2]=H[2,0]=-.25/130;H[1,2]=H[2,1]=-.32/10
    Ceng=np.linalg.solve(H,np.eye(6));CM=D@Ceng@D;HM=np.linalg.solve(D,np.linalg.solve(D,H).T).T
    rng=np.random.default_rng(20261103);e=rng.normal(size=6)*.01;e0=rng.normal(size=6)*.001
    metrics=dict(energy=rel((e-e0)@CM@(e-e0)-(D@(e-e0))@Ceng@(D@(e-e0)),(e-e0)@CM@(e-e0)),compliance=rel(CM@HM-np.eye(6),np.eye(6)))
    records.append(dict(case='Mandel_engineering_bridge',metrics=metrics))
    x=.31;z=.003;k=.4;J=np.array([[1-k*z,0,-k*x],[0,1,0],[k*x,0,1.]])
    w=np.array([1.,.2,.1]);A=axis(J,w);Q=rotation(.73);T=transform(A)
    metrics=dict(orthogonal=rel(A.T@A-np.eye(3),np.eye(3)),determinant=abs(float(np.linalg.det(A))-1),
        objective=rel(axis(Q@J,Q@w)-Q@A,A),Mandel_orthogonal=rel(T.T@T-np.eye(6),np.eye(6)),
        material_tensor=rel(T.T@mandel(tensor(e))-mandel(A.T@tensor(e)@A),e))
    records.append(dict(case='curved_reference_fiber_axis',metrics=metrics))
    rejected=False
    try:axis(np.eye(3),np.array([0.,0.,1.]))
    except ValueError:rejected=True
    records.append(dict(case='parallel_fiber_negative_control',metrics=dict(rejection_failure=float(not rejected))))
    ia=[0,1,3,4,5];ib=[2];Caa=CM[np.ix_(ia,ia)];Cab=CM[np.ix_(ia,ib)];Cba=Cab.T;Cbb=CM[np.ix_(ib,ib)]
    reduced=Caa-Cab@np.linalg.solve(Cbb,Cba);strain=e.copy();strain[ib]=e0[ib]-np.linalg.solve(Cbb,Cba@(e[ia]-e0[ia]));stress=CM@(strain-e0)
    by_compliance=np.linalg.solve(HM[np.ix_(ia,ia)],np.eye(5))
    ds=np.linalg.solve(HM[np.ix_(ia,ia)],e[ia]-e0[ia]);eb=e0[ib]+HM[np.ix_(ib,ia)]@ds
    metrics=dict(stress_normal=rel(stress[ib],stress),stress=rel(stress[ia]-reduced@(e[ia]-e0[ia]),stress),
        compliance_route=rel(reduced-by_compliance,reduced),recovery=rel(strain[ib]-eb,strain[ib]),
        energy=rel((strain-e0)@CM@(strain-e0)-(e[ia]-e0[ia])@reduced@(e[ia]-e0[ia]),(strain-e0)@CM@(strain-e0)))
    records.append(dict(case='normal_stress_Schur_and_recovery',metrics=metrics,omitting_Schur_detected=rel(Caa-reduced,reduced)>1e-4))
    layers=[(-.01,-.003,Ceng,rotation(.3),e0),(-.003,.01,.7*Ceng,rotation(-.61),-e0/2)]
    a=rng.normal(size=6)*.003;b=rng.normal(size=6)*.1;q=np.r_[a,b]
    K,load,offset=moments(layers);U,f=direct_stack(layers,a,b)
    records.append(dict(case='layer_moments_vs_tensor_quadrature',metrics=dict(energy=rel(U-(q@K@q/2-q@load+offset),U),force=rel(f-(K@q-load),f))))
    full=[(-.01,.01,Ceng,rotation(.3),e0)];split=[(-.01,-.002,Ceng,rotation(.3),e0),(-.002,.004,Ceng,rotation(.3),e0),(.004,.01,Ceng,rotation(.3),e0)]
    K0,l0,c0=moments(full);K1,l1,c1=moments(split);u0,f0=direct_stack(full,a,b);u1,f1=direct_stack(split,a,b)
    records.append(dict(case='homogeneous_layer_subdivision',metrics=dict(matrix=rel(K0-K1,K0),initial_load=rel(l0-l1,l0),constant=abs(c0-c1),
        energy=rel(u0-u1,u0),force=rel(f0-f1,f0))))
    for record in records:record['passed']=all(v<=1e-10 for v in record['metrics'].values()) and record.get('omitting_Schur_detected',True)
    sources=[Path(__file__),ROOT/'docs/specs/structure/material-layer-adoption-contract.md',ROOT/'rd/structure/material/records/material-layer-kernel-protocol.md',ROOT/'rd/structure/material/records/material-layer-kernel-preparation.json']
    dest.write_text(json.dumps(dict(records=records,all_checks_pass=all(r['passed'] for r in records),material_point_contract_only=True,
        element_plane_stress_adopted=False,real_material_validation=False,
        sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources}),indent=2),encoding='utf-8')
    print(json.dumps(records,indent=2))
if __name__=='__main__':main()
