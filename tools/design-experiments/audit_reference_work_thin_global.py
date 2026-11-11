"""Stored leading fields: tensor-energy, all shared orders, and supports."""
from pathlib import Path
import hashlib,json
import numpy as np
import reference_work_thin_limit_v3 as t
ROOT=Path(__file__).resolve().parents[2]
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def energy(E):return .3/(1.3*.4)*np.trace(E,axis1=1,axis2=2)**2+2/2.6*np.einsum('pij,pij->p',E,E)
def main():
    source=ROOT/'docs/verification/reference-work-thin-global-v4-result.json';parent=json.loads(source.read_text());dest=ROOT/'docs/verification/reference-work-thin-global-audit.json';assert not dest.exists();records=[];stops=[]
    for r in parent['records']:
        try:
            path=ROOT/r['artifact'];assert digest(path)==r['artifact_sha256'];data=np.load(path);total=np.zeros(2);reference=np.zeros(2);faces={};shared=support=constraint=0.
            for cell,(ids,coeff) in enumerate(zip(data['cells'],data['values'])):
                q=data['ref'][ids].reshape(-1);op=t.operators(q,15,9,True);s=op['state'];ns=s['nscalar'];v=coeff.reshape(3,ns,3,2).copy()
                for face in s['faces']:
                    Mc0,_=t.c.face_coefficients(q,face['axis'],face['sign'],False);T=face['Mcoeff'].T@Mc0;sl=slice(face['offset'],face['offset']+face['n'])
                    for order in range(3):v[order,sl]=np.einsum('ij,jak->iak',T,v[order,sl])
                for face,key in zip(s['faces'][2:],t.fit.edge_keys(cell,r['n'])):
                    sl=slice(face['offset'],face['offset']+face['n']);field=v[:,sl]
                    if key in faces:shared=max(shared,float(np.linalg.norm(field-faces[key])/max(1.,np.linalg.norm(field))))
                    else:faces[key]=field.copy()
                    if key[0]==0 and key[1]==0:support=max(support,float(np.linalg.norm(field)/max(1.,np.linalg.norm(v))))
                N,_=t.m.shape(op['p'],0,True);x=(N@q.reshape(75,3)[:25])[:,0];a=np.column_stack((np.ones(len(x)),np.zeros(len(x)),r['kappa']*x))
                for j,kind in enumerate(t.fit.KIND_TARGETS):
                    grad=[np.einsum('na,pjn->paj',v[k,:,:,j],s['Gwref']) for k in range(3)];C=op['C']
                    numerator=grad[0]@C[:,:,:,2]+grad[1]@C[:,:,:,1]+grad[2]@C[:,:,:,0];E=(numerator+numerator.transpose(0,2,1))/(2*op['d0'][:,None,None])
                    curvature=x*x if kind=='quartic' else np.exp(x);target=-op['p'][:,2,None,None]/2*curvature[:,None,None]*a[:,:,None]*a[:,None,:]/(1+r['kappa']**2*x*x)[:,None,None]
                    total[j]+=op['w']@(op['d0']*energy(E-target));reference[j]+=op['w']@(op['d0']*energy(target))
                    n0=grad[0]@C[:,:,:,0];n1=grad[0]@C[:,:,:,1]+grad[1]@C[:,:,:,0]
                    constraint=max(constraint,float(np.linalg.norm(n0+n0.transpose(0,2,1))/max(1.,np.linalg.norm(v))),float(np.linalg.norm(n1+n1.transpose(0,2,1))/max(1.,np.linalg.norm(v))))
                    cf,_=t.c.geometry_coefficients(np.asarray(q,dtype=np.float64).tobytes(),True)
                    for face in s['faces']:
                        j0=face['jump']@v[0,:,:,j];constraint=max(constraint,float(np.linalg.norm(j0)/max(1.,np.linalg.norm(v))))
                        j1=face['jump']@v[1,:,:,j]
                        if face['axis']==2:constraint=max(constraint,float(np.linalg.norm(j1)/max(1.,np.linalg.norm(v))))
                        else:
                            L,_=t.c.scalar_basis(face['points'],True);C1=np.einsum('pm,mab->pab',L,cf.reshape(len(t.c.IDS),3,3,5)[:,:,:,1],optimize=True);area=np.linalg.norm(C1[:,face['axis'],:],axis=1)
                            total[j]+=s['tau']*(face['w']*area)@np.sum(j1*j1,axis=1)
            errors=np.sqrt(total/reference);difference=max(abs(x-y) for x,y in zip(errors,r['primary']['error']))
            checks=dict(energy=difference<=1e-8,shared=shared<=1e-9,support=support<=1e-9,constraints=constraint<=1e-8)
            records.append(dict(kappa=r['kappa'],n=r['n'],amplitude=r['amplitude'],error=errors.tolist(),error_difference=difference,shared=shared,support=support,constraint=constraint,checks={k:bool(v) for k,v in checks.items()},passed=bool(all(checks.values()))))
        except Exception as exc:stops.append(dict(kappa=r['kappa'],n=r['n'],amplitude=r['amplitude'],reason=repr(exc)))
    paths=[Path(__file__),source,ROOT/'tools/design-experiments/reference_work_thin_limit_v3.py']
    dest.write_text(json.dumps(dict(records=records,stops=stops,all_checks_pass=len(records)==len(parent['records']) and not stops and all(r['passed'] for r in records),geometry_conditional=True,
        independent_global_resolve=False,sources={p.relative_to(ROOT).as_posix():digest(p) for p in paths}),indent=2),encoding='utf-8');print(json.dumps(dict(records=records,stops=stops),indent=2))
if __name__=='__main__':main()
