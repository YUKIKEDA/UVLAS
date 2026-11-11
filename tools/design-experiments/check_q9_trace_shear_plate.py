"""T-S36 global comparison after local, shared and independent gates; not product code."""
from pathlib import Path
import hashlib,json,math,platform,time
import numpy as np
import scipy
from scipy.sparse import coo_matrix
import check_q9_distorted_plate as base
import q9_trace_shear as projection
import check_q9_precision_diagnostic as precision
import check_q9_plate_convergence as plate
from check_q9_distorted_plate import mesh,patch,observables,digest,ROOT,F
PROTOCOL=F/'q9-trace-shear-plate-experiment.md'
def assemble(n,h,ref,nodes,cells,order,guard):
    values=[];f=np.zeros(len(nodes)*9)
    for _,_,ids,dofs in cells:
        guard();local=nodes[ids];values.append(projection.matrices(base.configuration(local),h,order)['Ktie'].reshape(-1))
        shape,x,w=base.surface_data(local,order)
        force=shape.T@(ref['pressure']*np.sin(np.pi*x[:,0])*np.sin(np.pi*x[:,1])*w)
        v=np.zeros((9,9));v[:,2]=force;v[:,8]=h*h/8*force;np.add.at(f,dofs,v.reshape(81))
    rows=np.concatenate([np.repeat(d,81) for _,_,_,d in cells]);cols=np.concatenate([np.tile(d,81) for _,_,_,d in cells])
    k=coo_matrix((np.concatenate(values),(rows,cols)),shape=(len(f),len(f))).tocsr();k.eliminate_zeros()
    return k,f

def reference(h):
    data=json.loads((F/'plate-modal-exact-reference-result.json').read_text());assert data['all_criteria_pass']
    return next(c['rounded_reference'] for c in data['cases'] if c['h']==h)
def run(n,h,out,stage_start):
    start=time.monotonic();deadline=min(start+240,stage_start+1200)
    def guard():
        if time.monotonic()>deadline:raise TimeoutError('case240_or_stage1200_seconds')
        if plate.memory_bytes()>2*1024**3:raise MemoryError('process_peak_2GiB')
        if sum(p.stat().st_size for p in out.rglob('*') if p.is_file())>64*1024**2:raise RuntimeError('output_64MiB')
    vertices,nodes,cells=mesh(n);geo=patch(nodes,cells,n)
    assert geo['polynomial_error']<=1e-11 and geo['area_error']<=1e-12 and geo['wrong_jacobian_error']>1e-3
    ref=reference(h);free=precision.free_dofs(n);matrices={};loads={};solves={};solutions={};records={}
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
    target=F/'q9-trace-shear-plate-result.json'
    if target.exists():raise FileExistsError(target)
    out=ROOT/'build/analysis'/('q9-trace-shear-plate-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir();start=time.monotonic();cases=[];checks=[];stop=None;uniform=False
    gates=['q9-trace-shear-local-result.json','q9-trace-shear-shared-result.json','q9-trace-shear-gates-audit.json','plate-modal-exact-reference-result.json']
    gate_hashes={n:digest(F/n) for n in gates}
    for n in gates:assert json.loads((F/n).read_text())['all_criteria_pass']
    audit=json.loads((F/gates[2]).read_text());assert audit['local_result_sha256']==gate_hashes[gates[0]] and audit['shared_result_sha256']==gate_hashes[gates[1]]
    for name in gates[:2]:
        data=json.loads((F/name).read_text());assert data['candidate_sha256']==digest(ROOT/'docs/q9-trace-shear-candidate.md')
        for n,h in data['dependency_sha256'].items():assert digest(ROOT/'tools/design-experiments'/n)==h
        for c in data['cases']:assert digest(ROOT/c['artifact'])==c['artifact_sha256']
    try:
        for h in [.01,.001]:
            for n in [2,4,8,16]:
                c=run(n,h,out,start);cases.append(c);print(json.dumps(dict(h=h,n=n,eligible=c['eligible'],observables=c['observables'],seconds=c['seconds'])),flush=True)
        for h in [.01,.001]:
            a=next(c for c in cases if c['h']==h and c['n']==8);b=next(c for c in cases if c['h']==h and c['n']==16)
            rate=math.log2(a['observables']['l2_relative']/b['observables']['l2_relative']) if a['observables'] and b['observables'] else None
            checks.append(dict(h=h,accuracy_pass=bool(b['eligible'] and all(b['observables'][k]<=.01 for k in ['center_relative','l2_relative','compliance_relative'])),last_l2_order=rate,order_pass=bool(a['eligible'] and b['eligible'] and (b['observables']['l2_relative']<=1e-8 or rate>=1.5))))
        a=next(c for c in cases if c['h']==.01 and c['n']==8);b=next(c for c in cases if c['h']==.001 and c['n']==8)
        uniform=bool(a['eligible'] and b['eligible'] and b['observables']['l2_relative']<=2*a['observables']['l2_relative']+1e-8)
    except Exception as exc:stop=repr(exc)
    result=dict(source_sha256=digest(__file__),protocol_sha256=digest(PROTOCOL),protocol=PROTOCOL.relative_to(ROOT).as_posix(),candidate_sha256=digest(ROOT/'docs/q9-trace-shear-candidate.md'),gate_sha256=gate_hashes,
        dependency_sha256={Path(m.__file__).name:digest(m.__file__) for m in [base,projection,precision,plate]},cases=cases,thickness_checks=checks,thin_uniformity_pass=uniform,stop_reason=stop,seconds=time.monotonic()-start,peak_bytes=plate.memory_bytes(),
        all_criteria_pass=bool(stop is None and len(cases)==8 and all(c['eligible'] for c in cases) and all(c['accuracy_pass'] and c['order_pass'] for c in checks) and uniform),product_shell_adopted=False)
    target.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(dict(all_criteria_pass=result['all_criteria_pass'],checks=checks,stop_reason=stop)))
    if not result['all_criteria_pass']:raise SystemExit(1)
if __name__=='__main__':main()
