"""Independent geometry, modal reference and observable audit of saved fields."""
from pathlib import Path
import hashlib,itertools,json,math,time
import numpy as np
from numpy.polynomial.legendre import leggauss
import check_q9_plate_convergence as memory
ROOT=Path(__file__).resolve().parents[2];F=ROOT/'docs/verification'
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
OFFSETS=[(-1,-1),(1,-1),(1,1),(-1,1),(0,-1),(1,0),(0,1),(-1,0),(0,0)]

def shape(x,y):
    def l(t,i):return .5*t*(t-1) if i==-1 else (1-t*t if i==0 else .5*t*(t+1))
    return np.array([l(x,i)*l(y,j) for i,j in OFFSETS])

def q4(x,y):
    w=np.array([(1-x)*(1-y),(1+x)*(1-y),(1+x)*(1+y),(1-x)*(1+y)])/4
    gx=np.array([-(1-y),1-y,1+y,-(1+y)])/4;gy=np.array([-(1-x),-(1+x),1+x,1-x])/4
    return w,np.column_stack((gx,gy))

def reference(h):
    p,w=leggauss(7);nu=.3;mu=1/(2*(1+nu));lam=nu/((1+nu)*(1-2*nu));k=np.zeros((4,4));pi=math.pi
    for t,weight in zip(p,w):
        z=t*h/2
        normal=np.array([[-pi*z,0,0,0],[0,-pi*z,0,0],[0,0,0,z]])
        shear=np.array([[0,1,pi,pi*z*z/2],[1,0,pi,pi*z*z/2],[pi*z,pi*z,0,0]])
        trace=normal.sum(axis=0)
        k+=(2*mu*(normal.T@normal)+lam*np.outer(trace,trace)+mu*shear.T@shear)*weight*h/8
    pressure=h**3/(12*(1-nu*nu))*(2*pi*pi)**2*(.001*h);rhs=pressure/4*np.array([0,0,1,h*h/8]);s=1/np.sqrt(k.diagonal())
    result=s*np.linalg.solve(s[:,None]*k*s[None,:],s*rhs)
    return dict(coefficients=result,pressure=pressure,center=float(result[2]),compliance=float(rhs@result))

def evaluate(n,vertices,nodes,q,f,ref):
    p,w=leggauss(11);error=0.;den=0.;center=None;minimum=math.inf;area=0.;node_error=0.
    for j in range(n):
        for i in range(n):
            corners=vertices[[j,j,j+1,j+1],[i,i+1,i+1,i]]
            ids=np.array([(2*j+y+1)*(2*n+1)+2*i+x+1 for x,y in OFFSETS])
            for index,(x,y) in zip(ids,OFFSETS):node_error=max(node_error,float(np.max(abs(nodes[index]-q4(x,y)[0]@corners))))
            for x,y in OFFSETS[:4]:minimum=min(minimum,4*n*n*float(np.linalg.det(corners.T@q4(x,y)[1])))
            for x,y in itertools.product(range(11),repeat=2):
                xi=p[x];eta=p[y];ww,gg=q4(xi,eta);position=ww@corners;measure=w[x]*w[y]*np.linalg.det(corners.T@gg)
                exact=ref['center']*math.sin(math.pi*position[0])*math.sin(math.pi*position[1]);actual=float(shape(xi,eta)@q[9*ids+2])
                error+=(actual-exact)**2*measure;den+=exact**2*measure;area+=measure
            # Independent Newton inverse; keep every valid cell on shared boundaries.
            parent=np.zeros(2)
            for _ in range(8):
                ww,gg=q4(*parent);parent-=np.linalg.solve(corners.T@gg,ww@corners-.5)
            if np.max(abs(parent))<=1+1e-10 and np.max(abs(q4(*parent)[0]@corners-.5))<=1e-12:
                value=float(shape(*parent)@q[9*ids+2])
                if center is not None:assert abs(value-center)<=1e-10*abs(ref['center'])
                center=value
    assert center is not None and node_error<=1e-13 and minimum>=.2 and abs(area-1)<=1e-12
    compliance=math.fsum(float(a)*float(b) for a,b in zip(f,q))
    return dict(center=center,compliance=compliance,center_relative=abs(center/ref['center']-1),compliance_relative=abs(compliance/ref['compliance']-1),l2_relative=math.sqrt(error/den))

def main():
    start=time.monotonic();source=F/'q9-distorted-plate-result.json';data=json.loads(source.read_text());observed=[];max_difference=0.;files=[]
    assert data['source_sha256']==digest(ROOT/'tools/design-experiments/check_q9_distorted_plate.py')
    assert data['protocol_sha256']==digest(ROOT/data['protocol']) and data['stop_reason'] is None
    assert data['controls']['all_criteria_pass'] and data['controls']['inverted_geometry_rejected']
    assert data['seconds']<=1200 and data['peak_bytes']<=2*1024**3
    for n,h in data['dependency_sha256'].items():assert digest(ROOT/'tools/design-experiments'/n)==h
    assert {(c['n'],c['h']) for c in data['cases']}=={(n,h) for n in [2,4,8,16] for h in [.01,.001]} and len(data['cases'])==8
    for c in data['cases']:
        assert c['seconds']<=240 and c['peak_bytes']<=2*1024**3 and c['eligible'] and c['assembly_pass']
        assert c['stiffness_7_vs_9']<=1e-8 and c['load_7_vs_9']<=1e-8 and c['symmetry']<=1e-12
        assert all(v<=1e-8 for d in c['precision_comparisons'].values() for v in d.values()) and all(v<=1e-5 for v in c['quadrature_response'].values())
        for group in c['solves'].values():
            for r in group.values():
                last=r['history'][-1];assert r['eligible'] and len(r['history'])<=6 and last['relative_residual']<=1e-9 and last['backward_normwise']<=1e-14 and last['backward_componentwise']<=1e-14
        file=ROOT/c['artifact'];assert digest(file)==c['artifact_sha256'];files.append(file)
        arrays=np.load(file);n=c['n'];h=c['h'];vertices=arrays['vertices'];nodes=arrays['nodes'];ref=reference(h)
        for j in range(n+1):
            for i in range(n+1):
                shift=.2/n*(-1)**(i+j)*math.sin(math.pi*i/n)*math.sin(math.pi*j/n) if 0<i<n and 0<j<n else 0
                assert np.max(abs(vertices[j,i]-[i/n+shift,j/n]))<=1e-13
        for k in ['pressure','center','compliance']:assert abs(ref[k]/c['reference'][k]-1)<=1e-10
        assert np.max(abs((ref['coefficients']-c['reference']['coefficients'])/ref['coefficients']))<=1e-10
        chosen=None
        for order in [7,9]:
            for digits in [70,90]:
                r=evaluate(n,vertices,nodes,arrays[f'q{order}_{digits}'],arrays[f'f{order}'],ref);saved=c['solves'][str(order)][str(digits)]['observables']
                for k in ['center','compliance']:assert abs((r[k]-saved[k])/ref[k])<=1e-10
                delta=abs(r['l2_relative']-saved['l2_relative']);max_difference=max(max_difference,delta)
                assert delta<=1e-8*max(saved['l2_relative'],1e-12)
                if order==9 and digits==90:chosen=r
        observed.append(dict(n=n,h=h,**chosen))
        assert time.monotonic()-start<=120 and memory.memory_bytes()<=512*1024**2
    checks=[]
    for h in [.01,.001]:
        a=next(c for c in observed if c['n']==8 and c['h']==h);b=next(c for c in observed if c['n']==16 and c['h']==h)
        rate=math.log2(a['l2_relative']/b['l2_relative']);old=next(c for c in data['thickness_checks'] if c['h']==h)
        accuracy=all(b[k]<=.01 for k in ['center_relative','l2_relative','compliance_relative']);order=b['l2_relative']<=1e-8 or rate>=1.5
        assert accuracy==old['accuracy_pass'] and order==old['order_pass'] and abs(rate-old['last_l2_order'])<=1e-8
        checks.append(dict(h=h,accuracy_pass=accuracy,order_pass=order,last_l2_order=rate))
    mid=next(c for c in observed if c['n']==8 and c['h']==.01);thin=next(c for c in observed if c['n']==8 and c['h']==.001)
    uniform=thin['l2_relative']<=2*mid['l2_relative']+1e-8;assert uniform==data['thin_uniformity_pass']
    acceptance=all(c['accuracy_pass'] and c['order_pass'] for c in checks) and uniform;assert acceptance==data['all_criteria_pass']
    assert sum(p.stat().st_size for p in files)<=64*1024**2
    result=dict(source_sha256=digest(__file__),memory_helper_sha256=digest(memory.__file__),audit_protocol_sha256=digest(F/'q9-distorted-plate-audit-protocol.md'),result_sha256=digest(source),
        independent_observables=observed,thickness_checks=checks,thin_uniformity_pass=uniform,thin_to_mid_l2_ratio=thin['l2_relative']/mid['l2_relative'],
        maximum_l2_difference=max_difference,seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),all_criteria_pass=True,experiment_acceptance_pass=acceptance,product_shell_adopted=False)
    (F/'q9-distorted-plate-artifact-audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ['all_criteria_pass','experiment_acceptance_pass','thin_to_mid_l2_ratio','maximum_l2_difference','seconds','peak_bytes']},indent=2))
if __name__=='__main__':main()
