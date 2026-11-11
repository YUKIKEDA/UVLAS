"""First bounded reference-only Ritz and physical-thickness checks."""
from pathlib import Path
import hashlib,json,time
import numpy as np
import pg_ritz_lift as lift
ROOT=Path(__file__).resolve().parents[2]

def main():
    dest=ROOT/'docs/verification/pg-ritz-reference-result.json';assert not dest.exists();records=[];stops=[]
    for n,parameters in [(2,(.4,0.,0.)),(2,(.4,.2,0.)),(2,(.4,-.2,.07)),(4,(.4,.2,0.))]:
        try:
            tick=time.monotonic();model=lift.prepare(n,parameters);print('PREPARED',n,parameters,time.monotonic()-tick,flush=True);p=np.array([[.13,.17],[.37,.31],[.61,.73],[.87,.91]]);v,d,ds=lift.evaluate(model,p);r,A,normal,dnormal,g=lift.geometry.graph(p,parameters);H=np.concatenate((A.transpose(0,2,1),normal[:,None,:]),axis=1);constraints=H@d+np.stack((np.einsum('pa,pan->pn',normal,v[1]),np.einsum('pa,pan->pn',normal,v[2]),np.zeros((len(p),32*n))),axis=1);derivative=0.;second=0.;step=1e-5
            for axis in [0,1]:
                samples=[]
                for factor in [-2,-1,1,2]:
                    point=p.copy();point[:,axis]+=factor*step;samples.append(lift.evaluate(model,point))
                fd=(samples[0][1]-8*samples[1][1]+8*samples[2][1]-samples[3][1])/(12*step);derivative=max(derivative,float(np.linalg.norm(fd-ds[axis])/max(1.,np.linalg.norm(ds[axis]))))
                for component,expected in [(1,3 if axis==0 else 4),(2,4 if axis==0 else 5)]:
                    fd=(samples[0][0][component]-8*samples[1][0][component]+8*samples[2][0][component]-samples[3][0][component])/(12*step);second=max(second,float(np.linalg.norm(fd-v[expected])/max(1.,np.linalg.norm(v[expected]))))
            thickness=0.;minimum=1.
            for h in [.01,.0001]:
                for variable in [False,True]:
                    top=lift.geometry.reference(p,np.ones(len(p)),h,parameters,variable);bottom=lift.geometry.reference(p,-np.ones(len(p)),h,parameters,variable);thickness=max(thickness,float(np.max(abs(np.linalg.norm(top[0]-bottom[0],axis=1)/top[3]-1))));minimum=min(minimum,float(top[2].min()),float(bottom[2].min()))
            checks=dict(stationarity=max(x['stationarity'] for x in model['metrics'])<=1e-8,condition=min(x['rcond'] for x in model['metrics'])>=1e-10,shear_normal=float(np.max(abs(constraints)))<=1e-9,second_derivative=second<=1e-8,director_derivative=derivative<=1e-8,normal_thickness=thickness<=1e-9,positive_reference=minimum>0)
            record=dict(n=n,parameters=parameters,metrics=model['metrics'],second_derivative=second,director_derivative=derivative,thickness=thickness,min_det=minimum,checks=checks,seconds=time.monotonic()-tick);records.append(record);print(record,flush=True)
            if not all(checks.values()):raise RuntimeError('Ritz_reference_gate')
        except Exception as exc:stops.append(dict(n=n,parameters=parameters,reason=repr(exc)));print('STOP',stops[-1],flush=True);break
    paths=[Path(__file__),Path(lift.__file__),Path(lift.geometry.__file__),ROOT/'docs/verification/pg-ritz-general-protocol.md',ROOT/'docs/verification/paired-geometry-ancf-physical-pde-audit.json'];result=dict(passed=len(records)==4 and not stops,records=records,stops=stops,R4_complete=False,main_shell_adopted=False,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths});dest.write_text(json.dumps(result,indent=2),encoding='utf-8')

if __name__=='__main__':main()
