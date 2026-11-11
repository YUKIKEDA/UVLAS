"""Separate bubble-gradient evaluation and trace/work checks for the flat candidate."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,itertools,json,time
import numpy as np
import q9_edge_trace_gradient as candidate
import audit_q9_curved_strain_screen as independent
import check_q9_curved_load as geometry
import check_q9_edge_moment_shared_space_v2 as shared
from q9_edge_moment_fields import quadrature
ROOT=rd_paths.ROOT;F=rd_paths.RECORDS
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    target=F/'q9-edge-trace-gradient-audit.json'
    if target.exists():raise FileExistsError(target)
    parent=F/'q9-edge-trace-gradient-gate-result.json';assert json.loads(parent.read_text())['all_criteria_pass']
    start=time.monotonic();cases=[]
    for name,h in itertools.product(['RECT','SKEW','TRAP'],[.1,.01]):
        q=geometry.geometry(name);state=candidate.prepare(q,h,0.,5);points,w=quadrature(7);J,E,B=candidate.evaluate(state,points);measure=w*np.linalg.det(J);V=measure.sum();psi=E[:,:,:81]@state['gradient_mode'];normal=state['gradient_mode'][72:75]
        _,_,grad=independent.data(q,h,points);gradient=grad[:,24,:];tensor=(normal[None,:,None]*gradient[:,None,:]+gradient[:,:,None]*normal[None,None,:])/2;other=np.einsum('kij,pij->pk',candidate.MB,tensor)
        independent_error=float(np.linalg.norm(psi-other)/np.linalg.norm(other));normpsi=np.sqrt(np.einsum('p,pi,pi->',measure,other,other));mean=float(np.linalg.norm(measure@other)/np.sqrt(V)/normpsi)
        normE=np.sqrt(np.einsum('p,pkn,pkn->n',measure,E,E));projection=np.einsum('p,pi,pin->n',measure,other,B-E);restore=float(np.max(abs(projection)/(normpsi*normE)))
        work=np.einsum('p,pkn->kn',measure,B-E);work_error=float(np.linalg.norm(work)/max(np.linalg.norm(np.einsum('p,pkn->kn',measure,E)),1e-30))
        boundary=[]
        for axis,sign in itertools.product(range(2),[-1,1]):
            ep=np.array([[x,sign,z] if axis==0 else [sign,x,z] for x,z in itertools.product([-.8,-.1,.6],[-.7,0.,.4])]);Je,Ee,_=candidate.evaluate(state,ep);pe=Ee[:,:,:81]@state['gradient_mode'];boundary.append(np.einsum('pi,kij,pk,pj->p',Je[:,:,axis],candidate.MB,pe,Je[:,:,2]*2/h))
        trace_error=float(np.linalg.norm(boundary)/max(np.linalg.norm(psi),1e-30));checks=dict(independent=independent_error<=1e-10,mean=mean<=1e-11,trace=trace_error<=1e-11,restore=restore<=1e-10,work7=work_error<=1e-10)
        cases.append(dict(geometry=name,h=h,independent=independent_error,mean=mean,trace=trace_error,restore=restore,work7=work_error,checks={k:bool(v) for k,v in checks.items()}))
    rejected=[]
    for name in ['SADDLE','NORMAL-C']:
        try:candidate.prepare(geometry.geometry(name),.1,0.,5)
        except ValueError:rejected.append(name)
    xy,cells,coords,_,_,_=shared.make_mesh('TRAP',2);states=[candidate.prepare(coords[ids].reshape(81),.01,0.,5) for _,_,ids,_ in cells];trace_errors=[]
    for j,i in itertools.product(range(2),repeat=2):
        index=2*j+i
        neighbors=([(index+1,0)] if i==0 else [])+([(index+2,1)] if j==0 else [])
        for neighbor,fixed_axis in neighbors:
            operators=[]
            for cell,sign in [(index,1),(neighbor,-1)]:
                ep=np.array([[sign,t,z] if fixed_axis==0 else [t,sign,z] for t,z in itertools.product([-.8,-.1,.6],[-.7,0.,.4])]);J,E,B=candidate.evaluate(states[cell],ep);row=2*np.einsum('pi,kij,pj->pk',J[:,:,1-fixed_axis],candidate.MB,J[:,:,2]*200);local=np.einsum('pk,pkn->pn',row,B[:,:,:81]);op=np.zeros((len(ep),9*len(xy)));dofs=(9*cells[cell][2][:,None]+np.arange(9)).reshape(-1);op[:,dofs]=local;operators.append(op)
            trace_errors.append(float(np.linalg.norm(operators[0]-operators[1])/max(np.linalg.norm(operators[0]),np.linalg.norm(operators[1]))))
    result=dict(source_sha256=digest(__file__),parent_sha256=digest(parent),protocol_sha256=digest(F/'q9-edge-trace-gradient-experiment.md'),dependencies={Path(m.__file__).name:digest(m.__file__) for m in [candidate,candidate.base,independent,geometry,shared]},cases=cases,curve_rejections=rejected,shared_trace_errors=trace_errors,seconds=time.monotonic()-start,peak_bytes=shared.memory.memory_bytes(),all_criteria_pass=bool(len(cases)==6 and all(all(c['checks'].values()) for c in cases) and len(rejected)==2 and len(trace_errors)==4 and max(trace_errors)<=1e-10),product_shell_adopted=False)
    target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps(result))
    if not result['all_criteria_pass']:raise SystemExit(1)
if __name__=='__main__':main()
