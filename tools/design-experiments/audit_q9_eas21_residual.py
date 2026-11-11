"""Independent strain/target, tensor Schur split, weighted thickness QR audit."""
from pathlib import Path
import hashlib,json,time
import numpy as np
from scipy.linalg import qr
import audit_q9_eas21 as source
import q9_developable_bending_approximation as mesh_source
ROOT=Path(__file__).resolve().parents[2]
MU=source.MU;LAM=source.LAM;LPS=2*MU*LAM/(LAM+2*MU)
CP=2*MU*np.eye(3)+LPS*np.outer([1.,1.,0.],[1.,1.,0.]);CH=np.linalg.cholesky(CP).T
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def quadratic(p):return 2*MU*np.sum(p*p,axis=-1)+LPS*(p[...,0]+p[...,1])**2
def local(error,weight,x,z,order):
    a=np.column_stack((np.ones_like(x),np.zeros_like(x),.4*x));a/=np.linalg.norm(a,axis=1)[:,None];n=np.cross(a,np.array([0.,1.,0.]));b=np.broadcast_to([0.,1.,0.],a.shape)
    def contraction(v,w):return np.einsum('pi,pij,pj->p',v,error,w)
    p=np.column_stack((contraction(a,a),contraction(b,b),np.sqrt(2.)*contraction(a,b)));normal=contraction(n,n)+LAM/(LAM+2*MU)*(p[:,0]+p[:,1]);shear=4*MU*(contraction(a,n)**2+contraction(b,n)**2)
    parts=np.array([0.,0.,0.,weight@((LAM+2*MU)*normal**2),weight@shear]);wm=weight.reshape(order,order,order);pm=p.reshape(order,order,order,3);zm=z.reshape(order,order,order);pzero=[];wzero=[]
    for i in range(order):
        for j in range(order):
            w=wm[i,j];X=np.column_stack((np.ones(order),zm[i,j]));Q,_=qr(X*np.sqrt(w)[:,None],mode='economic');Y=(pm[i,j]@CH.T)*np.sqrt(w)[:,None];coef=Q.T@Y;tail=Y-Q@coef;parts[:3]+=[coef[0]@coef[0],coef[1]@coef[1],np.sum(tail*tail)];pzero.append(w@pm[i,j]/w.sum());wzero.append(w.sum())
    pzero=np.array(pzero);wzero=np.array(wzero);mp=wzero@pzero/wzero.sum();memmean=float(wzero.sum()*quadratic(mp));mean=np.einsum('p,pij->ij',weight,error)/weight.sum();mean_energy=float(weight.sum()*source.energy(mean[None])[0]);total=float(weight@source.energy(error));return parts,mean_energy,memmean,total
def main():
    dest=ROOT/'docs/verification/q9-eas21-residual-audit.json'
    if dest.exists():raise FileExistsError(dest)
    start=time.monotonic();parent=ROOT/'docs/verification/q9-eas21-residual-result.json';p=json.loads(parent.read_text());shared=ROOT/'docs/verification/q9-eas21-shared-approximation-result.json';saved=json.loads(shared.read_text());records=[]
    for data in [p,saved]:
        for path,sha in data['sources'].items():assert digest(ROOT/path)==sha
    points,w=source.primary.quadrature(11);order=11
    for record,expected in zip(saved['records'],p['records']):
        path=ROOT/record['artifact'];assert digest(path)==record['artifact_sha256'];u=np.load(path)['u'];n,h=record['n'],record['h'];xy,cells,coords,_,_=mesh_source.mesh(.4,n,.2);parts=np.zeros(5);mean=memmean=total=ref=0.
        for _,_,ids,_ in cells:
            q=coords[ids].reshape(81);J,R=source.evaluate(source.prepare(q,h),points);phi,_,_=source.independent.data(q,h,points);z=h*points[:,2]/2;x=(phi@q.reshape(27,3))[:,0]/(1-.4*z);star=source.target(x,z);field=np.einsum('kij,pk->pij',source.MB,R@u[ids].reshape(81));weight=w*np.linalg.det(J);component,m,mm,t=local(field-star,weight,x,z,order);parts+=component;mean+=m;memmean+=mm;total+=t;ref+=weight@source.energy(star)
        expected=expected['routes']['saved_eas'];difference=float(np.max(abs(parts/ref-np.array(list(expected['energy_parts'].values())))));mean_difference=abs(mean/ref-expected['normalized_mean']);mem_difference=abs(memmean/ref-expected['normalized_membrane_mean']);error_difference=abs(np.sqrt(total/ref)-expected['error']);reconstruction=abs(parts.sum()-total)/total
        checks=dict(parts=difference<=1e-8,element_mean=mean_difference<=1e-8,membrane_mean=mem_difference<=1e-8,total=error_difference<=1e-8,reconstruction=reconstruction<=1e-10);assert all(checks.values());records.append(dict(n=n,h=h,parts_difference=difference,mean_difference=mean_difference,membrane_mean_difference=mem_difference,total_difference=error_difference,reconstruction=reconstruction,checks=checks));print(n,h,difference,error_difference,flush=True)
    sources=[Path(__file__),Path(source.__file__),Path(source.independent.__file__),parent,shared,ROOT/'docs/verification/q9-eas21-residual-experiment.md']
    dest.write_text(json.dumps(dict(records=records,all_checks_pass=True,independent_resolve=False,sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start),indent=2),encoding='utf-8')
if __name__=='__main__':main()
