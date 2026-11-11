"""Post-screen boundary values and matching two-element tensor traces."""
from pathlib import Path
import hashlib,itertools,json,time
import numpy as np
from numpy.polynomial.legendre import leggauss
import check_q9_curved_strain_screen as screen
import check_q9_curved_load as geometry
import check_q9_independent_load as common
import check_q9_plate_convergence as memory
ROOT=Path(__file__).resolve().parents[2];F=ROOT/'docs/verification'
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def maps(q,h,points,scheme):
    J,B=screen.kinematics(q,h,0.,points);ties,components=screen.trace.sampling_points();Js,Bs=screen.kinematics(q,h,0.,ties);Jc,_=screen.kinematics(q,h,0.,[[0,0,0]])
    A=J if scheme=='G-COV' else np.broadcast_to(Jc,J.shape);As=Js if scheme=='G-COV' else np.broadcast_to(Jc,Js.shape);S=screen.trace.interpolation(points)
    O=screen.voigt(screen.reduce(np.broadcast_to(screen.EB.transpose(1,2,0),(len(J),3,3,6)),np.broadcast_to(screen.EB.transpose(1,2,0),(36,3,3,6)),A,As,S,components))
    E=screen.reduce(screen.tensor(B),screen.tensor(Bs),A,As,S,components)
    reduced=np.einsum('pi,pijk,pj->pk',J[:,:,1],E,J[:,:,2]);raw=np.einsum('pi,pijk,pj->pk',J[:,:,1],screen.tensor(B),J[:,:,2])
    return J,O,reduced,raw
def two_elements(name):
    coords=[]
    for y in [-.5,0,.5]:
        for x in [-.5,-.25,0,.25,.5]:
            r=np.array([x,y,0.]);d=np.array([0.,0.,1.]);c=np.zeros(3)
            if name=='SKEW':r[0]+=.4*y
            elif name=='TRAP':r[0]*=1+.6*y
            elif name in ['SADDLE','NORMAL-C']:
                r[2]=.4*x*y;d=np.array([-.4*y,-.4*x,1.])
                if name=='NORMAL-C':d/=np.linalg.norm(d);c=.2*d
            coords.append([r,d,c])
    coords=np.array(coords);ids=[]
    for cell in [0,1]:ids.append(np.array([(int(t)+1)*5+2*cell+int(s)+1 for s,t in common.nodes]))
    return coords,ids
def main():
    target=F/'q9-curved-strain-boundary-result.json'
    if target.exists():raise FileExistsError(target)
    start=time.monotonic();parent=F/'q9-curved-strain-screen-result.json';assert __import__('json').loads(parent.read_text())['diagnostic_checks_pass'];records=[]
    p,_=leggauss(9);points=[]
    for axis,sign in itertools.product(range(2),[-1,1]):
        for t,z in itertools.product(p,[-1,0,1]):
            points.append([sign,t,z] if axis==0 else [t,sign,z])
    for name,h in itertools.product(['RECT','SKEW','TRAP','SADDLE','NORMAL-C'],[.1,.01]):
        q=geometry.geometry(name);coords,ids=two_elements(name);ql=coords[ids[0]].reshape(81);qr=coords[ids[1]].reshape(81);left=np.array([[1,t,z] for t,z in itertools.product([-.8,-.1,.6],[-.6,0,.7])]);right=left.copy();right[:,0]=-1
        for scheme in ['G-COV','C-COV']:
            _,O,_,_=maps(q,h,points,scheme);Jl,_,Ll,Rawl=maps(ql,h,left,scheme);Jr,_,Lr,Rawr=maps(qr,h,right,scheme)
            a=np.zeros((9,135));b=a.copy();ar=a.copy();br=a.copy();dl=(9*ids[0][:,None]+np.arange(9)).reshape(-1);dr=(9*ids[1][:,None]+np.arange(9)).reshape(-1);a[:,dl]=Ll;b[:,dr]=Lr;ar[:,dl]=Rawl;br[:,dr]=Rawr
            trace=float(np.linalg.norm(a-b)/max(np.linalg.norm(a),np.linalg.norm(b)));raw=float(np.linalg.norm(ar-br)/max(np.linalg.norm(ar),np.linalg.norm(br)));geodefect=float(np.max(abs(Jl[:,:,[1,2]]-Jr[:,:,[1,2]])))
            assert raw<=1e-10 and geodefect<=1e-12
            value=screen.field_error(O);records.append(dict(geometry=name,h=h,scheme=scheme,boundary_constant_value_error=value,shared_trace_error=trace,raw_trace_error=raw,geometry_difference=geodefect,value_gate=value<=1e-10,trace_gate=trace<=1e-10))
            assert time.monotonic()-start<=120 and memory.memory_bytes()<=512*1024**2
    result=dict(source_sha256=digest(__file__),parent_sha256=digest(parent),protocol_sha256=digest(F/'q9-curved-strain-boundary-experiment.md'),dependency_sha256={Path(m.__file__).name:digest(m.__file__) for m in [screen,geometry,common,memory]},cases=records,seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),diagnostic_checks_pass=True,operators_accepted=False,product_shell_adopted=False)
    target.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(cases=records,seconds=result['seconds'])))
if __name__=='__main__':main()
