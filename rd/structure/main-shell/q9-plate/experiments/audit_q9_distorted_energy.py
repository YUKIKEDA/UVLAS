"""Independent planar displacement gradients and energy integrals of saved fields."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,itertools,json,math,time
import numpy as np
from numpy.polynomial.legendre import leggauss
import check_q9_plate_convergence as memory
ROOT=rd_paths.ROOT;F=rd_paths.RECORDS
OFF=[(-1,-1),(1,-1),(1,1),(-1,1),(0,-1),(1,0),(0,1),(-1,0),(0,0)]
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def lag(x,i):
    return (.5*x*(x-1),x-.5) if i==-1 else ((1-x*x,-2*x) if i==0 else (.5*x*(x+1),x+.5))
def main():
    start=time.monotonic();data_path=F/'q9-distorted-energy-result.json';data=json.loads(data_path.read_text());parent_path=F/'q9-distorted-plate-result.json';parent=json.loads(parent_path.read_text())
    assert data['all_criteria_pass'] and data['parent_sha256']==digest(parent_path)
    assert data['source_sha256']==digest(ROOT/'rd/structure/main-shell/q9-plate/experiments/check_q9_distorted_energy.py')
    assert data['protocol_sha256']==digest(F/'q9-distorted-energy-experiment.md')
    assert data['parent_audit_sha256']==digest(F/'q9-distorted-plate-artifact-audit.json')
    for n,h in data['dependency_sha256'].items():assert digest(rd_paths.EXPERIMENTS/n)==h
    p,w=leggauss(11);t,wt=leggauss(3);N=[];ds=[];dt=[];weights=[];zz=[]
    for i,j,k in itertools.product(range(11),range(11),range(3)):
        a=[lag(p[i],s) for s,_ in OFF];b=[lag(p[j],s) for _,s in OFF]
        N.append([x[0]*y[0] for x,y in zip(a,b)]);ds.append([x[1]*y[0] for x,y in zip(a,b)]);dt.append([x[0]*y[1] for x,y in zip(a,b)])
        weights.append(w[i]*w[j]*wt[k]);zz.append(t[k])
    N=np.array(N);ds=np.array(ds);dt=np.array(dt);weights=np.array(weights);zz=np.array(zz);keys=['plane_even','plane_odd','thickness_residual','transverse_shear']
    mu=1/2.6;lam=.3/(1.3*.4)
    def integral(e,v):
        reflected=e.reshape(-1,3,6)[:,::-1,:].reshape(-1,6);even=(e+reflected)/2;odd=(e-reflected)/2
        def plane(x):return 2*mu*(x[:,0]**2+x[:,1]**2)+2*lam*mu/(lam+2*mu)*(x[:,0]+x[:,1])**2+mu*x[:,5]**2
        return np.array([v@plane(even),v@plane(odd),v@((lam+2*mu)*(e[:,2]+lam/(lam+2*mu)*(e[:,0]+e[:,1]))**2),v@(mu*(e[:,3]**2+e[:,4]**2))])
    maximum=0.;records=[]
    for c in parent['cases']:
        assert digest(ROOT/c['artifact'])==c['artifact_sha256'];arrays=np.load(ROOT/c['artifact']);n=c['n'];h=c['h'];z=h/2*zz;nodes=arrays['nodes'];q=arrays['q9_90'];A,B,C,H=c['reference']['coefficients']
        sums={name:np.zeros(4) for name in ['solution','reference','error','interpolation_error']}
        for j,i in itertools.product(range(n),repeat=2):
            ids=np.array([(2*j+b+1)*(2*n+1)+2*i+a+1 for a,b in OFF]);local=nodes[ids];xy=N@local
            jac=np.stack((ds@local,dt@local),axis=2);inverse=np.linalg.inv(jac);measure=weights*np.linalg.det(jac)*h/2
            sx=np.sin(np.pi*xy[:,0]);sy=np.sin(np.pi*xy[:,1]);cx=np.cos(np.pi*xy[:,0]);cy=np.cos(np.pi*xy[:,1]);s=sx*sy
            reference=np.column_stack((-np.pi*z*A*s,-np.pi*z*B*s,z*H*s,(B+np.pi*C+np.pi*H*z*z/2)*sx*cy,(A+np.pi*C+np.pi*H*z*z/2)*cx*sy,np.pi*z*(A+B)*cx*cy))
            qi=np.zeros((9,9));lx,ly=local.T;qi[:,2]=C*np.sin(np.pi*lx)*np.sin(np.pi*ly);qi[:,3]=A*np.cos(np.pi*lx)*np.sin(np.pi*ly);qi[:,4]=B*np.sin(np.pi*lx)*np.cos(np.pi*ly);qi[:,8]=H*np.sin(np.pi*lx)*np.sin(np.pi*ly)
            strains=[]
            for fields in [q.reshape(-1,9)[ids],qi]:
                r,d,c2=fields[:,:3],fields[:,3:6],fields[:,6:]
                us=ds@r+z[:,None]*(ds@d)+z[:,None]**2/2*(ds@c2);ut=dt@r+z[:,None]*(dt@d)+z[:,None]**2/2*(dt@c2)
                grad=np.einsum('pab,pbc->pac',np.stack((us,ut),axis=2),inverse);uz=N@d+z[:,None]*(N@c2)
                strains.append(np.column_stack((grad[:,0,0],grad[:,1,1],uz[:,2],uz[:,1]+grad[:,2,1],uz[:,0]+grad[:,2,0],grad[:,0,1]+grad[:,1,0])))
            for name,e in [('solution',strains[0]),('reference',reference),('error',strains[0]-reference),('interpolation_error',strains[1]-reference)]:sums[name]+=integral(e,measure)
        saved=next(r for r in data['cases'] if r['n']==n and r['h']==h);J=c['reference']['compliance'];relative=0.
        for name,values in sums.items():
            expected=np.array([saved['quadrature']['11']['components'][name][k] for k in keys]);relative=max(relative,float(np.max(abs(values-expected))/J))
        assert relative<=1e-8;maximum=max(maximum,relative);records.append(dict(n=n,h=h,maximum_component_difference_over_J=relative,components={name:dict(zip(keys,map(float,v))) for name,v in sums.items()}))
        assert time.monotonic()-start<=120 and memory.memory_bytes()<=512*1024**2
    result=dict(source_sha256=digest(__file__),result_sha256=digest(data_path),memory_helper_sha256=digest(memory.__file__),cases=records,
        maximum_component_difference_over_J=maximum,seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),all_criteria_pass=True,product_shell_adopted=False)
    (F/'q9-distorted-energy-artifact-audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ['all_criteria_pass','maximum_component_difference_over_J','seconds','peak_bytes']}))
if __name__=='__main__':main()
