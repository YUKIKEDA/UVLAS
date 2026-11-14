"""P-S36 global comparison after the fixed local gate; not product code."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,math,platform,time
import numpy as np
import scipy
from scipy.sparse import coo_matrix
import check_q9_distorted_plate as base
import check_q9_shear_projection_gate as projection
import check_q9_precision_diagnostic as precision
import check_q9_plate_convergence as plate
from check_q9_distorted_plate import mesh,patch,observables,digest,ROOT,F
PROTOCOL=F/'q9-shear-projection-plate-experiment.md'
def assemble(n,h,ref,nodes,cells,order,guard):
    values=[];f=np.zeros(len(nodes)*9)
    for _,_,ids,dofs in cells:
        guard();local=nodes[ids];values.append(projection.matrices(base.configuration(local),h,order)['Kp'].reshape(-1))
        shape,x,w=base.surface_data(local,order)
        force=shape.T@(ref['pressure']*np.sin(np.pi*x[:,0])*np.sin(np.pi*x[:,1])*w)
        v=np.zeros((9,9));v[:,2]=force;v[:,8]=h*h/8*force;np.add.at(f,dofs,v.reshape(81))
    rows=np.concatenate([np.repeat(d,81) for _,_,_,d in cells]);cols=np.concatenate([np.tile(d,81) for _,_,_,d in cells])
    k=coo_matrix((np.concatenate(values),(rows,cols)),shape=(len(f),len(f))).tocsr();k.eliminate_zeros()
    return k,f
def run(n,h,out,stage_start):
    start=time.monotonic();deadline=min(start+240,stage_start+1200)
    def guard():
        if time.monotonic()>deadline:raise TimeoutError('case240_or_stage1200_seconds')
        if plate.memory_bytes()>2*1024**3:raise MemoryError('process_peak_2GiB')
        if sum(p.stat().st_size for p in out.rglob('*') if p.is_file())>64*1024**2:raise RuntimeError('output_64MiB')
    vertices,nodes,cells=mesh(n);geo=patch(nodes,cells,n)
    assert geo['polynomial_error']<=1e-11 and geo['area_error']<=1e-12 and geo['wrong_jacobian_error']>1e-3
    ref=plate.reference(h);free=precision.free_dofs(n);matrices={};loads={};solves={};solutions={};records={}
    for order in [7,9]:
        k,f=assemble(n,h,ref,nodes,cells,order,guard);matrices[order]=k;loads[order]=f
    common_scale=1/np.sqrt(matrices[7].diagonal());from scipy.sparse import diags
    d=diags(common_scale);base=d@matrices[7]@d;difference=d@(matrices[9]-matrices[7])@d
    quadrature=float(np.linalg.norm(difference.data)/np.linalg.norm(base.data));load_difference=float(np.linalg.norm(loads[9]-loads[7])/np.linalg.norm(loads[9]))
    symmetry=max(float(np.linalg.norm((k-k.T).data)/np.linalg.norm(k.data)) for k in matrices.values())
    assembly_pass=quadrature<=1e-8 and load_difference<=1e-8 and symmetry<=1e-12 and all(np.all(k.diagonal()>0) for k in matrices.values())
    if assembly_pass:
        for order in [7,9]:
            guard();a,b,scale,lu=precision.scaled_system(matrices[order],loads[order],free)
            solves[order]={}
            for digits in [70,90]:
                guard();r,q=precision.multiprecision(a,b,scale,lu,free,len(loads[order]),digits,deadline)
                r['observables']=observables(n,vertices,nodes,cells,q,loads[order],ref)
                solves[order][digits]=r;solutions[order,digits]=q
            records[str(order)]=precision.difference(solutions[order,70],solutions[order,90],solves[order][70]['observables'],solves[order][90]['observables'])
        integration_response=precision.difference(solutions[7,90],solutions[9,90],solves[7][90]['observables'],solves[9][90]['observables'])
        eligible=all(r['eligible'] for group in solves.values() for r in group.values()) and all(v<=1e-8 for group in records.values() for v in group.values()) and all(v<=1e-5 for v in integration_response.values())
    else:eligible=False;integration_response=None
    file=out/f'n{n}-h{h:g}.npz'
    np.savez_compressed(file,nodes=nodes,vertices=vertices,f7=loads[7],f9=loads[9],**{f'q{o}_{p}':q for (o,p),q in solutions.items()})
    guard()
    return dict(n=n,h=h,reference=ref,geometry=geo,stiffness_7_vs_9=quadrature,load_7_vs_9=load_difference,symmetry=symmetry,assembly_pass=bool(assembly_pass),solves=solves,
        precision_comparisons=records,quadrature_response=integration_response,eligible=bool(eligible),observables=solves[9][90]['observables'] if assembly_pass else None,
        seconds=time.monotonic()-start,peak_bytes=plate.memory_bytes(),artifact=file.relative_to(ROOT).as_posix(),artifact_sha256=digest(file))

def main():
    target=F/'q9-shear-projection-plate-result.json'
    if target.exists():raise FileExistsError(target)
    out=ROOT/'build/analysis'/('q9-shear-projection-plate-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir()
    start=time.monotonic();cases=[];stopped=None;control=None;checks=[]
    try:
        gate=F/'q9-shear-projection-gate-result.json';control=json.loads(gate.read_text());assert control['all_criteria_pass']
        assert control['source_sha256']==digest(projection.__file__) and control['candidate_sha256']==digest(ROOT/'rd/structure/main-shell/q9-plate/q9-shear-projection-candidate.md')
        for record in control['cases']:assert digest(ROOT/record['artifact'])==record['artifact_sha256']
        for h in [.01,.001]:
            for n in [2,4,8,16]:
                result=run(n,h,out,start);cases.append(result)
                print(json.dumps(dict(h=h,n=n,eligible=result['eligible'],observables=result['observables'],seconds=result['seconds'])),flush=True)
        for h in [.01,.001]:
            a=next(r for r in cases if r['h']==h and r['n']==8);b=next(r for r in cases if r['h']==h and r['n']==16)
            rate=math.log2(a['observables']['l2_relative']/b['observables']['l2_relative']) if a['observables'] and b['observables'] else None
            checks.append(dict(h=h,accuracy_pass=bool(b['eligible'] and all(b['observables'][k]<=.01 for k in ['center_relative','l2_relative','compliance_relative'])),last_l2_order=rate,
                order_pass=bool(a['eligible'] and b['eligible'] and (b['observables']['l2_relative']<=1e-8 or rate>=1.5))))
        mid=next(r for r in cases if r['h']==.01 and r['n']==8);thin=next(r for r in cases if r['h']==.001 and r['n']==8)
        uniform=bool(mid['eligible'] and thin['eligible'] and thin['observables']['l2_relative']<=2*mid['observables']['l2_relative']+1e-8)
    except Exception as exc:stopped=repr(exc);uniform=False
    data=dict(protocol=PROTOCOL.relative_to(ROOT).as_posix(),protocol_sha256=digest(PROTOCOL),source_sha256=digest(__file__),
        dependency_sha256={Path(m.__file__).name:digest(m.__file__) for m in [base,projection,plate,precision]},python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
        controls=control,gate_sha256=digest(F/'q9-shear-projection-gate-result.json'),baseline_sha256=digest(F/'q9-distorted-plate-result.json'),candidate_sha256=digest(ROOT/'rd/structure/main-shell/q9-plate/q9-shear-projection-candidate.md'),cases=cases,thickness_checks=checks,thin_uniformity_pass=uniform,stop_reason=stopped,seconds=time.monotonic()-start,peak_bytes=plate.memory_bytes(),
        all_criteria_pass=bool(stopped is None and len(cases)==8 and all(r['eligible'] for r in cases) and all(r['accuracy_pass'] and r['order_pass'] for r in checks) and uniform),
        product_shell_adopted=False,limitation='P-S36 projected transverse shear; flat isotropic linear comparison only, not general curved/anisotropic/nonlinear adoption')
    target.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8');print(json.dumps(dict(all_criteria_pass=data['all_criteria_pass'],checks=checks,stop_reason=stopped)))
    if not data['all_criteria_pass']:raise SystemExit(1)
if __name__=='__main__':main()
