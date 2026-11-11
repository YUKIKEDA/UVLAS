"""Freeze eight geometry-bound dual maps before implementing the new solve probe."""
from pathlib import Path
import hashlib,json,math,time
import numpy as np
import avl_symmetric_observation_weights as weights
import avl_symmetry_reconstruction as expansion
import check_avl_conservative_observation as obs
import check_q9_plate_convergence as memory
ROOT=Path(__file__).resolve().parents[2];F=ROOT/'docs/verification'
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    output=F/'avl-symmetric-solve-weight-result.json'
    if output.exists():raise FileExistsError(output)
    begin=time.monotonic();folder=ROOT/'build/analysis/avl-symmetric-solve-weights-v2';folder.mkdir();unit=weights.manufactured_geometry()
    target_path=F/'avl-tip-redistribution-result.json';target=json.loads(target_path.read_text());assert target['methods_criteria_pass']
    small_path=F/'avl-span1280-capacity-ks-runs.json';small=json.loads(small_path.read_text());assert small['completed']
    records=[];selections=[]
    for case in ['FW-V01','FW-V02']:
        run=next(r for r in small['runs'] if r['case']==case and r['m']==6)
        geometry=folder/(case+'-small-geometry');expansion.expand(ROOT/run['directory']/'out',6,40,geometry)
        _,bands,_=obs.native_model(dict(case=case,m=6,n=40,directory=geometry.relative_to(ROOT).as_posix()),'small-weight-check')
        selections.append(dict(case=case,allocation='small',m=6,n=40,raw_directory=run['directory']+'/out',derived_directory=geometry.relative_to(ROOT).as_posix(),bands=obs.public(bands)))
    selections+=target['meshes']
    for rec in selections:
        geometry=ROOT/rec['derived_directory'];rows=np.loadtxt(geometry/'uv-rows.txt',skiprows=2);forces=np.loadtxt(geometry/'uv-forces.txt',skiprows=2)
        finite=forces[forces[:,2]==0];assert np.array_equal(finite[:,0],np.arange(1,len(rows)+1));points=finite[:,3:6]
        full,reduced,order=weights.weights(rows,points,rec['case'],rec['m']);native=len(rows)//2
        q0=rows[native:,13]/math.sin(math.pi/180);expected=weights.public_vector(rec['bands']);actual=weights.apply(reduced,q0)
        original_error=float(np.max(abs(actual-expected)/np.maximum(1.,abs(expected))));assert original_error<=1e-12
        trials=[]
        for label,q in [('constant',np.ones(native)),('variable',np.sin(np.arange(native)*.37)+.1*np.cos(np.arange(native)*.11))]:
            g=np.concatenate([q[order],q]);mapped=weights.apply(reduced,q);direct=weights.direct(rows,points,g,rec['case'],rec['m'])
            error=float(np.max(abs(mapped-direct)/np.maximum(1.,abs(direct))));fold=float(np.max(abs(mapped-weights.apply(full,g))/np.maximum(1.,abs(mapped))))
            symmetry=0.
            for k in range(5):
                v=mapped[k*20:(k+1)*20];symmetry=max(symmetry,float(np.max(abs(v-((-1 if k==1 else 1)*v[::-1])))))
            dual=np.cos(np.arange(100)*.3)/np.maximum(1.,np.linalg.norm(reduced,axis=1))
            lhs=math.fsum(dual*mapped);rhs=math.fsum((reduced.T@dual)*q);transpose=abs(lhs-rhs)/max(1.,abs(lhs))
            assert max(error,fold,symmetry,transpose)<=1e-12
            trials.append(dict(kind=label,direct_error=error,fold_error=fold,symmetry_error=symmetry,transpose_error=transpose))
        archive=folder/f'{rec["case"]}-{rec["allocation"]}-weights.npz';widths=np.array([b['width'] for b in rec['bands']])
        np.savez_compressed(archive,matrix=reduced,q0=q0,widths=widths,original_observations=expected)
        record=dict(case=rec['case'],allocation=rec['allocation'],m=rec['m'],n=rec['n'],native_unknowns=native,raw_directory=rec['raw_directory'],
            geometry_directory=rec['derived_directory'],geometry_files={n:digest(geometry/n) for n in ['uv-rows.txt','uv-forces.txt']},
            input=rec['raw_directory'].rsplit('/out',1)[0]+'/wing.avl',weight_path=archive.relative_to(ROOT).as_posix(),weight_sha256=digest(archive),
            original_reproduction_error=original_error,trials=trials,all_criteria_pass=True)
        records.append(record);print(json.dumps(dict(case=rec['case'],allocation=rec['allocation'],weights=[100,native],seconds=time.monotonic()-begin)),flush=True)
        if time.monotonic()-begin>600 or memory.memory_bytes()>512*1024**2:raise RuntimeError('weight_resource_limit')
    total=sum(p.stat().st_size for p in folder.rglob('*') if p.is_file());assert total<32*1024**2
    result=dict(protocol='docs/verification/avl-symmetric-solve-sensitivity-experiment.md',protocol_sha256=digest(F/'avl-symmetric-solve-sensitivity-experiment.md'),
        source_sha256=digest(__file__),dependency_sha256={Path(m.__file__).name:digest(m.__file__) for m in [weights,expansion,obs,memory]},
        parent_sha256={target_path.name:digest(target_path),small_path.name:digest(small_path)},manufacturing=unit,records=records,
        layout=weights.LAYOUT,seconds=time.monotonic()-begin,peak_bytes=memory.memory_bytes(),output_bytes=total,all_criteria_pass=len(records)==8 and unit['all_criteria_pass'])
    obs.write(output,result);print(json.dumps(dict(all_criteria_pass=result['all_criteria_pass'],records=len(records),bytes=total)))
if __name__=='__main__':main()
