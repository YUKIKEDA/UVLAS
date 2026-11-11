"""Stored full fields: independent work, position error, traces and sensitivity."""
from pathlib import Path
import hashlib,json,time
import numpy as np
from scipy.linalg import cholesky,solve_triangular,qr
import mixed_mean_stress_thin_v6 as t
ROOT=Path(__file__).resolve().parents[2]
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    source=ROOT/'docs/verification/mixed-mean-stress-thin-v6-result.json';parent=json.loads(source.read_text());dest=ROOT/'docs/verification/mixed-mean-stress-thin-audit.json';assert not dest.exists();records=[];stops=[];start=time.monotonic()
    for record in parent['records']:
        try:
            path=ROOT/record['artifact'];assert digest(path)==record['artifact_sha256'];data=np.load(path);shared={};share=support=work=jump=fluct=0.;err=np.zeros(3);ref=np.zeros(3);sensitivity=[]
            for index,(ids,v) in enumerate(zip(data['cells'],data['values'])):
                q=data['reference'][ids].reshape(-1);op=t.old.operators(q,15,9,True);s=op['state'];ns=s['nscalar'];nc=3*ns;fields=v.reshape(2,ns,3,3).copy()
                for face in s['faces']:
                    Mc0,_=t.c.face_coefficients(q,face['axis'],face['sign'],False);T=face['Mcoeff'].T@Mc0;sl=slice(face['offset'],face['offset']+face['n']);fields[:,sl]=np.einsum('ij,ojak->oiak',T,fields[:,sl])
                for face,key in zip(s['faces'][2:],t.fit.edge_keys(index,record['n'])):
                    sl=slice(face['offset'],face['offset']+face['n']);value=fields[:,sl]
                    if key in shared:share=max(share,float(np.linalg.norm(value-shared[key])/max(1.,np.linalg.norm(value))))
                    else:shared[key]=value.copy()
                    if key[0]==0 and key[1]==0:support=max(support,float(np.linalg.norm(value)/max(1.,np.linalg.norm(fields))))
                for face in s['faces']:
                    jump=max(jump,float(np.linalg.norm(np.einsum('pn,nak->pak',face['jump'],fields[0]))/max(1.,np.linalg.norm(fields))))
                    if face['axis']==2:jump=max(jump,float(np.linalg.norm(np.einsum('pn,nak->pak',face['jump'],fields[1]))/max(1.,np.linalg.norm(fields))))
                E=(op['A'][1]@fields[0].reshape(nc,3)+op['A'][0]@fields[1].reshape(nc,3))/op['d0'][:,None,None];plane,wp=t.m.quadrature(15,2);_,wz=np.polynomial.legendre.leggauss(9);E=E.reshape(len(plane),9,6,3);mean=np.einsum('z,pzkt->pkt',wz/2,E);fluct=max(fluct,float(np.linalg.norm(E-mean[:,None])/max(1.,np.linalg.norm(E))))
                L,meta=t.stress.evaluate(q,record['kappa'],plane,True);weight=2*wp*op['d0'].reshape(len(plane),9)[:,0];gram=np.einsum('p,pir,ij,pjs->rs',weight,L,np.linalg.inv(t.m.C),L,optimize=True);moment=np.einsum('p,pkr,pkt->rt',weight,L,mean,optimize=True);res=solve_triangular(cholesky(gram,lower=True),moment,lower=True);work=max(work,float(np.linalg.norm(res)/max(1.,np.linalg.norm(fields))))
                Q,R,coef=t.scalar_position(s);I=Q@R;a=np.linalg.lstsq(I,fields[0].reshape(ns,-1),rcond=1e-12)[0];p=np.column_stack((plane,np.zeros(len(plane))));Lc,_=t.c.scalar_basis(p,True);position=(Lc@coef@a).reshape(len(p),3,3);target=t.target(q,p,record['kappa']);N,g=t.m.shape(p,0,True);J=np.einsum('na,pnb->pab',q.reshape(75,3)[:25],g);area=np.linalg.norm(np.cross(J[:,:,0],J[:,:,1]),axis=1);err+=np.einsum('p,pak,pak->k',wp*area,position-target,position-target);ref+=np.einsum('p,pak,pak->k',wp*area,target,target)
                if record['local'][index]['q01'].get('block_equilibrated'):
                    primary=t.local(q,record['kappa']);ind=t.local(q,record['kappa'],15,9,True);P=qr(primary['Z'],mode='economic')[0];QI=qr(ind['Z'],mode='economic')[0];sensitivity.append(dict(cell=index,primary_rank=primary['info']['ranks'],independent_rank=ind['info']['ranks'],span_difference=float(np.linalg.norm(P-QI@(QI.T@P))),independent_basis_digits=meta['exact_space_certificate']['basis_preprocessing_decimal_digits']))
            errors=np.sqrt(err/ref);difference=max(abs(errors-np.array(record['error'])));checks=dict(position_error=difference<=1e-8,work=work<=1e-9,shared=share<=1e-9,support=support<=1e-9,stabilization=jump<=1e-9,thickness_fluctuation=fluct<=1e-9)
            records.append(dict(kappa=record['kappa'],n=record['n'],amplitude=record['amplitude'],error=errors.tolist(),difference=float(difference),work=work,shared=share,support=support,stabilization=jump,fluctuation=fluct,completion_space_sensitivity=sensitivity,checks=checks,passed=all(checks.values())));print('AUDIT',record['kappa'],record['n'],record['amplitude'],checks,'sensitivity',sensitivity,flush=True)
        except Exception as exc:stops.append(dict(kappa=record['kappa'],n=record['n'],amplitude=record['amplitude'],reason=repr(exc)));print('AUDIT STOP',repr(exc),flush=True)
    paths=[Path(__file__),source,ROOT/'tools/design-experiments/mixed_mean_stress_thin_v6.py',ROOT/'tools/design-experiments/mixed_mean_stress_space_v4.py',ROOT/'tools/design-experiments/mixed_mean_stress_basis.py']
    dest.write_text(json.dumps(dict(records=records,stops=stops,all_audit_checks=len(records)==len(parent['records']) and not stops and all(r['passed'] for r in records),independent_global_resolve=False,seconds=time.monotonic()-start,sources={p.relative_to(ROOT).as_posix():digest(p) for p in paths}),indent=2),encoding='utf-8')
if __name__=='__main__':main()
