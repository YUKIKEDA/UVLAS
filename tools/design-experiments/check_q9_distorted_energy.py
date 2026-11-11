"""Positive energy partition of frozen solutions; no repeated equilibrium solve."""
from pathlib import Path
import hashlib,json,math,time
import numpy as np
import check_q9_distorted_plate as experiment
import check_q9_independent_load as material
import check_q9_plate_convergence as memory
ROOT=Path(__file__).resolve().parents[2];F=ROOT/'docs/verification'
KEYS=['plane_even','plane_odd','thickness_residual','transverse_shear']
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def reference_strain(x,coefficients):
    A,B,C,H=coefficients;xx,yy,z=x.T;sx=np.sin(np.pi*xx);sy=np.sin(np.pi*yy);cx=np.cos(np.pi*xx);cy=np.cos(np.pi*yy)
    return np.column_stack((-np.pi*z*A*sx*sy,-np.pi*z*B*sx*sy,z*H*sx*sy,
        (B+np.pi*C+.5*np.pi*H*z*z)*sx*cy,(A+np.pi*C+.5*np.pi*H*z*z)*cx*sy,np.pi*z*(A+B)*cx*cy))

def strain(grad,q):
    g=np.einsum('ia,pib->pab',q.reshape(27,3),grad)
    return np.column_stack((g[:,0,0],g[:,1,1],g[:,2,2],g[:,1,2]+g[:,2,1],g[:,0,2]+g[:,2,0],g[:,0,1]+g[:,1,0]))

def partition(e,measure):
    lam=material.lam;mu=material.mu
    def plane(v):
        return 2*mu*(v[:,0]**2+v[:,1]**2)+2*lam*mu/(lam+2*mu)*(v[:,0]+v[:,1])**2+mu*v[:,5]**2
    reflected=e.reshape(-1,3,6)[:,::-1,:].reshape(-1,6)
    even=(e+reflected)/2;odd=(e-reflected)/2
    values=np.array([measure@plane(even),measure@plane(odd),measure@((lam+2*mu)*(e[:,2]+lam/(lam+2*mu)*(e[:,0]+e[:,1]))**2),measure@(mu*(e[:,3]**2+e[:,4]**2))])
    direct=float(measure@np.einsum('pi,ij,pj->p',e,material.elastic,e))
    return values,direct,float(measure@(plane(e)-plane(even)-plane(odd)))

def interpolant(nodes,coefficients):
    A,B,C,H=coefficients;x,y=nodes.T;sx=np.sin(np.pi*x);sy=np.sin(np.pi*y);q=np.zeros((len(nodes),9))
    q[:,2]=C*sx*sy;q[:,3]=A*np.cos(np.pi*x)*sy;q[:,4]=B*sx*np.cos(np.pi*y);q[:,8]=H*sx*sy
    return q.reshape(-1)

def integrate(c,arrays,order,guard):
    nodes=arrays['nodes'];q=arrays['q9_90'];coeff=c['reference']['coefficients'];qi=interpolant(nodes,coeff)
    totals={name:np.zeros(4) for name in ['solution','reference','error','interpolation_error']};direct={name:0. for name in totals};cross={name:0. for name in totals}
    scalar,g,weights=experiment.table(order,c['h'])
    for _,_,ids,dofs in experiment.plate.connectivity(c['n']):
        guard();geometry=experiment.configuration(nodes[ids]).reshape(27,3)
        jac=np.einsum('ia,pib->pab',geometry,g);grad=np.linalg.solve(jac.transpose(0,2,1),g.transpose(0,2,1)).transpose(0,2,1)
        measure=weights*np.linalg.det(jac);x=scalar@geometry;actual=strain(grad,q[dofs]);ref=reference_strain(x,coeff);interp=strain(grad,qi[dofs])
        for name,e in [('solution',actual),('reference',ref),('error',actual-ref),('interpolation_error',interp-ref)]:
            value,full,odd=partition(e,measure);totals[name]+=value;direct[name]+=full;cross[name]+=odd
    return dict(components={name:dict(zip(KEYS,map(float,values))) for name,values in totals.items()},direct=direct,odd_cross=cross)

def main():
    output=F/'q9-distorted-energy-result.json'
    if output.exists():raise FileExistsError(output)
    parent_path=F/'q9-distorted-plate-result.json';parent=json.loads(parent_path.read_text());audit_path=F/'q9-distorted-plate-artifact-audit.json';audit=json.loads(audit_path.read_text())
    assert audit['all_criteria_pass'] and audit['result_sha256']==digest(parent_path)
    assert parent['source_sha256']==digest(experiment.__file__)
    for n,h in parent['dependency_sha256'].items():assert digest(ROOT/'tools/design-experiments'/n)==h
    start=time.monotonic();records=[];stop=None
    def guard():
        if time.monotonic()-start>120:raise TimeoutError('120second_diagnostic')
        if memory.memory_bytes()>512*1024**2:raise MemoryError('512MiB_diagnostic')
    try:
        for c in parent['cases']:
            assert c['eligible'] and digest(ROOT/c['artifact'])==c['artifact_sha256']
            arrays=np.load(ROOT/c['artifact']);runs={o:integrate(c,arrays,o,guard) for o in [9,11]};chosen=runs[11];J=c['reference']['compliance']
            work=math.fsum(float(f)*float(q) for f,q in zip(arrays['f9'],arrays['q9_90']))
            summation=max(abs(sum(v.values())-chosen['direct'][name])/J for name,v in chosen['components'].items())
            quadrature=max(abs(runs[9]['components'][name][k]-chosen['components'][name][k])/J for name in chosen['components'] for k in KEYS)
            cross=max(abs(v)/J for v in chosen['odd_cross'].values());reference_identity=abs(chosen['direct']['reference']/J-1)
            equilibrium_identity=abs(chosen['direct']['solution']-work)/J;galerkin_identity=abs(chosen['direct']['error']-(J-work))/J
            checks=dict(positive=all(v>=0 for values in chosen['components'].values() for v in values.values()),reconstruction=summation<=1e-8,
                integration=quadrature<=1e-7,odd_even=cross<=1e-8,reference=reference_identity<=1e-8,equilibrium=equilibrium_identity<=1e-7,galerkin=galerkin_identity<=1e-7)
            normalized={name:{k:v/J for k,v in values.items()} for name,values in chosen['components'].items()}
            shares={name:{k:v/sum(values.values()) for k,v in values.items()} for name,values in chosen['components'].items()}
            records.append(dict(n=c['n'],h=c['h'],artifact_sha256=c['artifact_sha256'],quadrature=runs,normalized=normalized,shares=shares,
                reference_J=J,work=work,defects=dict(reconstruction=summation,integration=quadrature,odd_even=cross,reference=reference_identity,equilibrium=equilibrium_identity,galerkin=galerkin_identity),checks=checks))
            print(json.dumps(dict(n=c['n'],h=c['h'],checks=checks,error_shares=shares['error'],normalized_error=sum(normalized['error'].values()),normalized_interpolation_error=sum(normalized['interpolation_error'].values()))),flush=True)
    except Exception as exc:stop=repr(exc)
    result=dict(protocol_sha256=digest(F/'q9-distorted-energy-experiment.md'),source_sha256=digest(__file__),parent_sha256=digest(parent_path),parent_audit_sha256=digest(audit_path),
        dependency_sha256={Path(m.__file__).name:digest(m.__file__) for m in [experiment,material,memory]},cases=records,stop_reason=stop,seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),
        all_criteria_pass=bool(stop is None and len(records)==8 and all(all(r['checks'].values()) for r in records)),product_shell_adopted=False)
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');assert output.stat().st_size<=8*1024**2
    if not result['all_criteria_pass']:raise SystemExit(1)
if __name__=='__main__':main()
