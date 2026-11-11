"""Vandermonde plus physical-z complex-step check of shared thin constraints."""
from pathlib import Path
import hashlib,json,itertools,time
import numpy as np
import q3_curved_independent as independent
import q3_thin_expansion as primary
ROOT=Path(__file__).resolve().parents[2]
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def expansion(q,points):
    N,g=independent.shape(points);z=1e-25j;phi_g=np.zeros((len(points),48,3),dtype=complex)
    for node in range(16):
        for power,factor in [(0,1.),(1,1.),(2,.5)]:
            col=node*3+power;phi_g[:,col,:2]=g[:,node]*factor*z**power
            if power:phi_g[:,col,2]=N[:,node]*factor*power*z**(power-1)
    J=np.einsum('na,pnb->pab',q.reshape(48,3),phi_g);grad=phi_g@np.linalg.inv(J);B=np.einsum('kij,pnj->pkni',primary.element.MB,grad).reshape(-1,6,144);return J.real,B.real,B.imag/1e-25
def main():
    dest=ROOT/'docs/verification/q3-thin-limit-audit.json'
    if dest.exists():raise FileExistsError(dest)
    start=time.monotonic();parent=ROOT/'docs/verification/q3-thin-limit-result.json';data=json.loads(parent.read_text());assert data['all_checks_pass'];records=[];p,w=np.polynomial.legendre.leggauss(11);points=np.array(list(itertools.product(p,p)));weights=np.array([a*b for a,b in itertools.product(w,w)]);N,g=independent.shape(points)
    for path,sha in data['sources'].items():assert digest(ROOT/path)==sha
    for record in data['records']:
        path=ROOT/record['artifact'];assert digest(path)==record['artifact_sha256'];saved=np.load(path);free=saved['free'];lookup={int(v):i for i,v in enumerate(free)};maxB0=maxB1=0.;kernel_square=0.;matrix_square=0.;out={kind:dict(error=0.,reference=0.,zero=0.,h1=0.) for kind in ['quadratic','quartic']};area=0.
        for ids in saved['ids']:
            q=saved['coords'][ids].reshape(144);J,B0,B1=expansion(q,points);_,P0,P1=primary.expansion(q,points);maxB0=max(maxB0,float(np.linalg.norm(B0-P0)/np.linalg.norm(P0)));maxB1=max(maxB1,float(np.linalg.norm(B1-P1)/np.linalg.norm(P1)));weight=weights*np.linalg.det(J);nodeXY=N@saved['xy'][ids];x=nodeXY[:,0];D=1+record['kappa']**2*x*x;a=np.column_stack((np.ones(len(x)),np.zeros(len(x)),record['kappa']*x));kernel_local=np.zeros((144,saved['kernel'].shape[1]))
            for local,node in enumerate(ids):
                if int(node) in lookup:kernel_local[9*local:9*local+6]=saved['kernel'][6*lookup[int(node)]:6*lookup[int(node)]+6]
            residual=B0@kernel_local;kernel_square+=float(np.sum(weight[:,None,None]*residual*residual));matrix_square+=float(np.sum(weight[:,None,None]*B0[:,:,primary.RD]**2));Jsurf=np.einsum('na,pnb->pab',saved['coords'][ids,0],g);metric=Jsurf.transpose(0,2,1)@Jsurf;aw=weights*np.sqrt(np.linalg.det(metric));area+=aw.sum()
            for kind in out:
                u=saved['u_'+kind][ids];v=np.ones_like(x) if kind=='quadratic' else x*x;target=-v[:,None,None]*a[:,:,None]*a[:,None,:]/D[:,None,None];field=np.einsum('kij,pk->pij',primary.element.MB,B1@u.reshape(144));e=field-target;lam=.3/(1.3*.4);mu=1/2.6
                def energy(E):return lam*np.trace(E,axis1=1,axis2=2)**2+2*mu*np.einsum('pij,pij->p',E,E)
                out[kind]['error']+=weight@energy(e);out[kind]['reference']+=weight@energy(target);out[kind]['zero']+=float(np.sum(weight[:,None]*(B0@u.reshape(144))**2));val=N@u[:,3:9];du=np.einsum('na,pnb->pab',u[:,3:9],g);out[kind]['h1']+=aw@(np.sum(val*val,axis=1)+np.einsum('pai,pij,paj->p',du,np.linalg.inv(metric),du))
        kernel_residual=float(np.sqrt(kernel_square/matrix_square)/np.linalg.norm(saved['kernel']));observations={}
        for kind,value in out.items():
            error=float(np.sqrt(value['error']/value['reference']));h1=float(np.sqrt(value['h1']/area));difference=abs(error-record['solutions'][kind]['observed_error']);h1difference=abs(h1/record['solutions'][kind]['observed_h1']-1);zero=float(np.sqrt(value['zero']/value['reference']));checks=dict(error=difference<=1e-8,h1=h1difference<=1e-8,zero=zero<=1e-9);assert all(checks.values());observations[kind]=dict(error=error,error_difference=difference,h1_difference=h1difference,zero_constraint=zero,checks=checks)
        proportional=None
        if record['kappa']==.4 and record['amplitude']==.2:
            reference=primary.field('quadratic',.4,saved['xy'])[free,:6].reshape(-1);K=saved['kernel'];assert K.shape[1]==1;coef=(K[:,0]@reference)/(K[:,0]@K[:,0]);proportional=float(np.linalg.norm(K[:,0]*coef-reference)/np.linalg.norm(reference));assert proportional<=1e-8
        assert maxB0<=1e-10 and maxB1<=1e-10 and kernel_residual<=1e-9;records.append(dict(kappa=record['kappa'],amplitude=record['amplitude'],n=record['n'],operator0_difference=maxB0,operator1_difference=maxB1,kernel_residual=kernel_residual,quadratic_mode_difference=proportional,observations=observations));print(record['kappa'],record['amplitude'],record['n'],'B',maxB0,maxB1,'mode',proportional,flush=True)
    sources=[Path(__file__),Path(independent.__file__),Path(primary.__file__),parent,ROOT/'docs/verification/q3-thin-limit-experiment.md']
    dest.write_text(json.dumps(dict(records=records,all_checks_pass=True,independent_global_resolve=False,finite_thickness_lower_bound=False,sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start),indent=2),encoding='utf-8')
if __name__=='__main__':main()
