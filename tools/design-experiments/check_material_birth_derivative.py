from pathlib import Path
import json,hashlib,numpy as np
from numpy.polynomial.legendre import leggauss
import uvlm_wake_assembly_v2 as aero
import uvlm_surface_bindings as binding
ROOT=Path(__file__).resolve().parents[2]
def main():
    b=binding.q1('wing',1,1,[[0,-.5,0],[0,.5,0],[1,.5,0],[1,-.5,0]]);model=aero.Model([b],.025);flow=aero.Flow([2,0,.12]);q=b['q0'];g=model.geometry(q);info=model.surface_info[0];te=g['X'][info['te']];X=te[None]+np.arange(3)[:,None,None]*.15*flow.U;rows=(('wing','attached'),('wing','old',1),('wing','old',2));ids=(tuple(model.vertex_ids[i] for i in info['te']),)+tuple(tuple((row,j) for j in range(2)) for row in rows[1:]);w=aero.Wake('wing',X,np.zeros((3,1)),ids,rows,(0.,-.15,-.3));p=g['P'][0];n=g['n'][0];B,W,_,_=aero.field_basis(model,g,(w,),p);coeff=np.array([1.,.8,.6]);gamma=np.linalg.solve((n@(B+W@coeff[:,None]@model.T))[None,:],np.array([-n@flow.U]));w=aero.Wake('wing',X,(coeff*gamma[0])[:,None],ids,rows,(0.,-.15,-.3));state=aero.State(model.identity,1,0.,q,np.zeros_like(q),gamma,(w,),gamma,.15,'explicit_state',flow.identity,'free',None);V=aero.convection(model,state,flow,'free')[0];base=n@(aero.field(model,g,(w,),gamma,p)[0]+flow.U);assert abs(base)<1e-12
    z,weights=leggauss(24);flux=np.zeros(3);span=te[1]-te[0]
    for zz,ww in zip((z+1)/2,weights/2):
        pos=(1-zz)*te[0]+zz*te[1];vel=(1-zz)*V[0,0]+zz*V[0,1];area=np.cross(vel,span);R=p-pos;rr=R@R;D=rr+model.eps**2;flux+=ww*gamma[0]*((rr-2*model.eps**2)*area-3*R*(R@area))/(4*np.pi*D**2.5)
    old=aero.field(model,g,(w,),gamma,p,dwakes=[V])[1];gd=np.array([.17]);acc=np.tile([0,0,.05],(model.ns,1));expected=float(n@(B@gd+old+flux-model.Nc[0]@acc));wrong=float(n@((B+W[:,:1]@model.T)@gd+old-model.Nc[0]@acc));errors=[]
    def residual(h):
        qq=q+.5*h*h*acc;gg=model.geometry(qq,h*acc);newX=np.concatenate([gg['X'][info['te']][None],X+h*V]);release=('wing','new',h);newids=(ids[0],tuple((release,j) for j in range(2)))+ids[1:];newrows=(rows[0],release)+rows[1:];gamma_new=gamma+h*gd;mu=np.concatenate([(model.T@gamma_new).reshape(1,1),w.mu]);wake=aero.Wake('wing',newX,mu,newids,newrows,(h,0.,-.15,-.3));return float(gg['n'][0]@(aero.field(model,gg,(wake,),gamma_new,gg['P'][0])[0]+flow.value(gg['P'][0],h)-gg['CV'][0]))
    for h in [.002,.001,.0005,.00025]:
        derivative=(4*residual(h/2)-residual(h)-3*base)/h;errors.append(abs(derivative-expected))
    order=float(np.log2(errors[-2]/errors[-1]));assert errors[-1]<1e-6*max(abs(expected),1) and 1.7<order<2.3 and abs(wrong-expected)>1e-3;paths=[Path(__file__),ROOT/'docs/j1-initialization-output-derivative-addendum.md',ROOT/'docs/verification/uvlm-wake-adoption-closure.json'];result=dict(passed=True,scope='actual weak right derivative of newborn strip; old material coefficients preserved',analytic_derivative=expected,finite_right_derivative_errors=errors,final_order=order,wrong_Aeff_and_omitted_flux_error=abs(wrong-expected),sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths});dest=ROOT/'docs/verification/material-birth-derivative-result.json';assert not dest.exists();dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
