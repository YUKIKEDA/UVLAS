"""Independent tensor energy/work and all-order traces of stored bending fields."""
from pathlib import Path
import hashlib,json,time
import numpy as np
from scipy.linalg import cholesky,solve_triangular
import mixed_mean_stress_bending as t
ROOT=Path(__file__).resolve().parents[2]
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def tensor_energy(E):return .3/(1.3*.4)*np.trace(E,axis1=0,axis2=1)**2+2/2.6*np.einsum('ij,ij->',E,E)
def energy(E):return .3/(1.3*.4)*np.trace(E,axis1=1,axis2=2)**2+2/2.6*np.einsum('pij,pij->p',E,E)
def phi(p):
    a=np.polynomial.legendre.legvander(p[:,0],5);b=np.polynomial.legendre.legvander(p[:,1],5);return np.column_stack([a[:,i]*b[:,j]*np.sqrt((2*i+1)*(2*j+1)) for i in range(6) for j in range(6)])
def main():
    source=ROOT/'docs/verification/mixed-mean-stress-bending-result.json';parent=json.loads(source.read_text());dest=ROOT/'docs/verification/mixed-mean-stress-bending-audit.json';assert not dest.exists();records=[];stops=[];start=time.monotonic()
    for r in parent['records']:
        try:
            path=ROOT/r['artifact'];assert digest(path)==r['artifact_sha256'];data=np.load(path);total=np.zeros(3);reference=np.zeros(3);errpos=np.zeros(3);refpos=np.zeros(3);faces={};shared=support=work=zero=jump=position_reconstruction=0.
            for cell,(ids,coeff) in enumerate(zip(data['cells'],data['values'])):
                q=data['ref'][ids].reshape(-1);op=t.operators(q,15,9,True);s=op['state'];ns=s['nscalar'];v=coeff.reshape(3,ns,3,3).copy()
                for face in s['faces']:
                    M0,_=t.c.face_coefficients(q,face['axis'],face['sign'],False);T=face['Mcoeff'].T@M0;sl=slice(face['offset'],face['offset']+face['n']);v[:,sl]=np.einsum('ij,ojak->oiak',T,v[:,sl])
                norm=max(1.,np.linalg.norm(v[:2]))
                for face,key in zip(s['faces'][2:],t.fit.edge_keys(cell,r['n'])):
                    sl=slice(face['offset'],face['offset']+face['n']);field=v[:,sl]
                    if key in faces:shared=max(shared,float(np.linalg.norm(field-faces[key])/max(1.,np.linalg.norm(field))))
                    else:faces[key]=field.copy()
                    if key[0]==0 and key[1]==0:support=max(support,float(np.linalg.norm(field)/max(1.,np.linalg.norm(v))))
                plane,wp=t.m.quadrature(15,2);pp=np.column_stack((plane,np.zeros(len(plane))));corners=q.reshape(75,3)[[0,4,20,24]];xp=t.mixed.stress.p.bilinear(corners[:,0]);yp=t.mixed.stress.p.bilinear(corners[:,1])
                def ev(poly,p):return sum(float(value)*p[:,0]**i*p[:,1]**j for (i,j),value in poly.items()) if poly else np.zeros(len(p))
                x=ev(xp,op['p']);xx=ev(xp,pp);jac=ev(t.mixed.stress.p.deriv(xp,0),pp)*ev(t.mixed.stress.p.deriv(yp,1),pp)-ev(t.mixed.stress.p.deriv(xp,1),pp)*ev(t.mixed.stress.p.deriv(yp,0),pp);area=jac*np.sqrt(1+r['kappa']**2*xx*xx)
                I=np.zeros((ns,36));I[:76]=s['V'].T@((s['w']/8)[:,None]*phi(s['p']))
                for face in s['faces']:I[face['offset']:face['offset']+face['n']]=face['hat'].T@((face['w']/4)[:,None]*phi(face['points']))
                aa=np.linalg.lstsq(I,v[0].reshape(ns,-1),rcond=1e-12)[0];position_reconstruction=max(position_reconstruction,float(np.linalg.norm(I@aa-v[0].reshape(ns,-1))/norm));position=(phi(pp)@aa).reshape(len(pp),3,3);btarget=np.column_stack((xx**4/12,xx**5,np.exp(xx)-1-xx));ctarget=np.column_stack((-r['kappa']*xx**5/15,-5*r['kappa']*xx**6/6,-r['kappa']*((xx-1)*np.exp(xx)+1-xx*xx/2)));ptarget=np.zeros_like(position);ptarget[:,0]=ctarget;ptarget[:,2]=btarget;errpos+=np.einsum('p,pak,pak->k',wp*area,position-ptarget,position-ptarget);refpos+=np.einsum('p,pak,pak->k',wp*area,ptarget,ptarget)
                L,_=t.mixed.stress.evaluate(q,r['kappa'],plane,True);ww=2*wp*op['d0'].reshape(len(plane),9)[:,0];gram=np.einsum('p,pir,ij,pjs->rs',ww,L,np.linalg.inv(t.m.C),L,optimize=True);chol=cholesky(gram,lower=True)
                af=np.column_stack((np.ones(len(x)),np.zeros(len(x)),r['kappa']*x));cf,_=t.c.geometry_coefficients(np.asarray(q,dtype=np.float64).tobytes(),True)
                for j,curvature in enumerate([x*x,20*x**3,np.exp(x)]):
                    grad=[np.einsum('na,pjn->paj',v[k,:,:,j],s['Gwref']) for k in range(3)];C=op['C'];n2=grad[0]@C[:,:,:,2]+grad[1]@C[:,:,:,1]+grad[2]@C[:,:,:,0];E=(n2+n2.transpose(0,2,1))/(2*op['d0'][:,None,None]);target=-op['p'][:,2,None,None]/2*curvature[:,None,None]*af[:,:,None]*af[:,None,:]/(1+r['kappa']**2*x*x)[:,None,None];total[j]+=op['w']@(op['d0']*energy(E-target));reference[j]+=op['w']@(op['d0']*energy(target))
                    n0=grad[0]@C[:,:,:,0];zero=max(zero,float(np.linalg.norm(n0+n0.transpose(0,2,1))/norm));n1=grad[0]@C[:,:,:,1]+grad[1]@C[:,:,:,0];E0=(n1+n1.transpose(0,2,1))/(2*op['d0'][:,None,None]);mean=np.einsum('z,pzij->pij',np.polynomial.legendre.leggauss(9)[1]/2,E0.reshape(len(plane),9,3,3));em=np.einsum('kij,pij->pk',t.m.MB,mean);moment=np.einsum('p,pkr,pk->r',ww,L,em);work=max(work,float(np.linalg.norm(solve_triangular(chol,moment,lower=True))/norm))
                    for face in s['faces']:
                        j0=face['jump']@v[0,:,:,j];jump=max(jump,float(np.linalg.norm(j0)/norm));j1=face['jump']@v[1,:,:,j]
                        if face['axis']==2:jump=max(jump,float(np.linalg.norm(j1)/norm))
                        else:
                            Lf,_=t.c.scalar_basis(face['points'],True);C1=np.einsum('pm,mab->pab',Lf,cf.reshape(len(t.c.IDS),3,3,5)[:,:,:,1],optimize=True);measure=np.linalg.norm(C1[:,face['axis']],axis=1);total[j]+=s['tau']*(face['w']*measure)@np.sum(j1*j1,axis=1)
            errors=np.sqrt(total/reference);positions=np.sqrt(errpos/refpos);de=float(np.max(abs(errors-np.array(r['primary']['error']))));dp=float(np.max(abs(positions-np.array(r['primary']['position_error']))));checks={key:bool(value) for key,value in dict(energy=de<=1e-8,position=dp<=1e-8,shared=shared<=1e-9,support=support<=1e-9,work=work<=1e-9,zero=zero<=1e-9,stabilization=jump<=1e-9,position_reconstruction=position_reconstruction<=1e-9).items()};records.append(dict(kappa=r['kappa'],n=r['n'],amplitude=r['amplitude'],error=errors.tolist(),position_error=positions.tolist(),energy_difference=de,position_difference=dp,shared=shared,support=support,work=work,zero=zero,jump=jump,position_reconstruction=position_reconstruction,checks=checks,passed=all(checks.values())));print('AUDIT',r['kappa'],r['n'],r['amplitude'],checks,flush=True)
        except Exception as exc:stops.append(dict(kappa=r['kappa'],n=r['n'],amplitude=r['amplitude'],reason=repr(exc)));print('AUDIT STOP',repr(exc),flush=True)
    paths=[Path(__file__),source,ROOT/'tools/design-experiments/mixed_mean_stress_bending.py',ROOT/'docs/verification/mixed-mean-stress-bending-protocol.md']
    dest.write_text(json.dumps(dict(records=records,stops=stops,all_audit_checks=parent['all_numerical_checks'] and len(records)==5 and not stops and all(r['passed'] for r in records),independent_global_resolve=False,geometry_conditional=True,seconds=time.monotonic()-start,sources={p.relative_to(ROOT).as_posix():digest(p) for p in paths}),indent=2),encoding='utf-8')
if __name__=='__main__':main()
