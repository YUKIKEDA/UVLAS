import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import numpy as np
import hashlib,json,time,traceback,sys
import structural_layer_damping as law
import connection_shell as pg
import material_equivalent_input as material
ROOT=rd_paths.ROOT
def rel(a,b):return float(np.linalg.norm(a)/max(np.linalg.norm(b),1e-30))
def main():
    start=time.monotonic();rows=[];failures=[]
    try:
        model=pg.shell.prepare(1,(.4,.2,.07),.1,False,'extended');model=dict(model,physical_scale=1.,physical_h=.1,span=1.,width=1.);d=pg.data(model,order=29);q0=d['q0'];rng=np.random.default_rng(1515);u=q0@np.array([[.03,.01,0],[0,-.02,.01],[.02,0,.04]]).T;v=q0@np.array([[.01,-.02,0],[.02,.01,.01],[0,-.01,-.02]]).T;theta=.2+.4*d['p'][:,0];axes=np.zeros((len(theta),3,3));axes[:,0,0]=axes[:,1,1]=np.cos(theta);axes[:,1,0]=np.sin(theta);axes[:,0,1]=-np.sin(theta);axes[:,2,2]=1.;T=material.transform(axes);lower=d['p'][:,2]<0;diagonal=np.where(lower[:,None],np.array([.3,.4,.7,.2,.8,.5]),np.array([.8,.6,.4,.5,.1,.3]));factor=np.sqrt(diagonal)[:,:,None]*np.swapaxes(T,1,2);f,Kq,Kv,info=law.shell_damping(d,u,v,factor);direction=q0@np.array([[.02,.01,-.02],[-.01,.02,.01],[.01,.02,.01]]).T;du=law.shell_damping(d,u+1e-25j*direction,v,factor,False)[0].imag/1e-25;dv=law.shell_damping(d,u,v+1e-25j*direction,factor,False)[0].imag/1e-25;qerr=rel(du.ravel()-Kq@direction.ravel(),du);verr=rel(dv.ravel()-Kv@direction.ravel(),dv);physical=axes.transpose(0,2,1)@((info['F'].transpose(0,2,1)@info['V']+info['V'].transpose(0,2,1)@info['F'])/2)@axes;rate=np.einsum('kab,pab->pk',material.MB,physical);independent=float(np.sum(d['weights'][:,None]*diagonal*rate*rate));power=abs(independent-np.sum(f*v))/independent;assert max(qerr,verr,power)<1e-8;rows.append(dict(case='anisotropic_layerwise_viscosity_in_reference_material_axes',state_tangent_error=qerr,velocity_tangent_error=verr,independent_material_power_error=float(power),two_distinct_layers=True,varying_reference_axes=True))
        total=np.zeros_like(f)
        for mask in [lower,~lower]:
            part={key:np.asarray(d[key])[mask] for key in ['G','weights']};total+=law.shell_damping(part,u,v,factor[mask],False)[0]
        add=rel(total-f,f);assert add<1e-10;rows.append(dict(case='layer_partition_additivity',force_error=add))
    except Exception as exc:failures.append(dict(reason=repr(exc),traceback=traceback.format_exc()));print('FAIL',failures[-1],flush=True)
    print(json.dumps(rows,indent=2),flush=True)
    if '--publish' in sys.argv:
        dest=ROOT/'rd/structure/integration/records/structural-layer-damping-result.json';assert not dest.exists();paths=[Path(__file__),Path(law.__file__),Path(material.__file__),ROOT/'rd/structure/integration/experiments/structural_load_damping.py',ROOT/'rd/structure/integration/records/structural-final-integration-protocol.md'];dest.write_text(json.dumps(dict(passed=len(rows)==2 and not failures,records=rows,failures=failures,seconds=time.monotonic()-start,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}),indent=2),encoding='utf-8');print('STRUCTURAL_LAYER_DAMPING_RESULT',len(rows)==2 and not failures,flush=True)
if __name__=='__main__':main()
