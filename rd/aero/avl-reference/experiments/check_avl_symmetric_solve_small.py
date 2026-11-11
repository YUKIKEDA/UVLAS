"""Small native-half matrix, transpose, Decimal80 residual and observation gate."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
from decimal import Decimal,localcontext
import hashlib,json,math,time
import numpy as np
import avl_symmetric_observation_weights as weights
ROOT=rd_paths.ROOT;F=rd_paths.RECORDS
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    output=F/'avl-symmetric-solve-small-check.json'
    if output.exists():raise FileExistsError(output)
    begin=time.monotonic();run_path=F/'avl-symmetric-solve-small-runs.json';raw=json.loads(run_path.read_text());assert raw['completed'] and len(raw['runs'])==2
    records=[]
    for r in raw['runs']:
        folder=ROOT/r['directory'];n=r['native_unknowns'];assert n==120
        assert digest(folder/'uq/uq-matrix.bin')==digest(folder/'out/uh-matrix.bin')
        a=np.fromfile(folder/'uq/uq-matrix.bin',dtype='<f8').reshape((n,n),order='F');v=np.loadtxt(folder/'uq/uq-vectors.txt',skiprows=1)
        assert v.shape==(n,10) and np.array_equal(v[:,0],np.arange(1,n+1))
        rhs=v[:,1];q=v[:,2:6];res=v[:,6:10];native=np.loadtxt(folder/'out/uv-rows.txt',skiprows=2)
        assert np.array_equal(rhs,-native[:,6])
        direct=np.linalg.solve(a,rhs);direct_error=float(np.max(abs(direct-q[:,1]))/max(1.,float(np.max(abs(direct)))))
        recurrence_error=max(float(np.max(abs(q[:,0]+np.linalg.solve(a,res[:,0])-q[:,2]))),float(np.max(abs(q[:,2]+np.linalg.solve(a,res[:,2])-q[:,3]))))
        trans=np.loadtxt(folder/'uq/uq-transpose.txt',skiprows=1);transpose_error=0.
        for k in [1,2,3]:
            part=trans[trans[:,1]==k];assert part.shape==(n,5);expected=np.linalg.solve(a.T,part[:,2])
            transpose_error=max(transpose_error,float(np.max(abs(expected-part[:,3]))/max(1.,float(np.max(abs(expected))))))
        cross=np.loadtxt(folder/'uq/uq-residual-cross.txt',skiprows=1);assert cross.shape==(4*n,5)
        quad_error=0.;comp_error=0.;scale_error=0.
        with localcontext() as ctx:
            ctx.prec=80;ad=[[Decimal.from_float(float(x)) for x in row] for row in a];qd=[[Decimal.from_float(float(x)) for x in col] for col in q.T]
            for row in cross:
                i=int(row[0])-1;k=int(row[1])-1;b=Decimal.from_float(float(rhs[i]));products=[x*y for x,y in zip(ad[i],qd[k])]
                exact=b-sum(products,Decimal(0));scale=abs(b)+sum((abs(x) for x in products),Decimal(0))
                quad_error=max(quad_error,abs(float(exact)-row[2]));comp_error=max(comp_error,abs(float(exact)-row[3])/max(1.,float(scale)))
                scale_error=max(scale_error,abs(float(scale)-row[4])/max(1.,float(scale)));assert row[3]==res[i,k]
        cache=np.load(ROOT/r['weight_path']);original=weights.apply(cache['matrix'],q[:,0]);observation_error=float(np.max(abs(original-cache['original_observations'])/np.maximum(1.,abs(cache['original_observations']))))
        assert direct_error<=1e-10 and recurrence_error<=1e-10 and transpose_error<=1e-10 and max(quad_error,comp_error,scale_error,observation_error)<=1e-12
        assert all(r['parent_outputs_equal'].values()) and all(r['before_after_equal'].values()) and r['probe_meta']['STATE_UNCHANGED']==1
        records.append(dict(case=r['case'],matrix_bytes_equal=True,direct_error=direct_error,recurrence_error=recurrence_error,transpose_error=transpose_error,
            decimal_quad_absolute_difference=quad_error,decimal_compensated_scaled_difference=comp_error,scale_error=scale_error,observation_error=observation_error,
            numpy_condition_2=float(np.linalg.cond(a)),all_criteria_pass=True))
    result=dict(source_sha256=digest(__file__),weights_module_sha256=digest(weights.__file__),small_runs_sha256=digest(run_path),
        preparation_sha256=digest(F/'avl-symmetric-solve-preparation.json'),records=records,seconds=time.monotonic()-begin,all_criteria_pass=len(records)==2 and all(r['all_criteria_pass'] for r in records),
        reference_acceptance=False)
    output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
