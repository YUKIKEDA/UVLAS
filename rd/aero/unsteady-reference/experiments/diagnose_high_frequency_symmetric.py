import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import sys,json,hashlib,math,numpy as np
import run_high_frequency_final_bounded as wrapper
from check_vortex_core_load import kernel
ROOT=rd_paths.ROOT
def unique_field(g,X,gamma,P,N,eps):
    m,n=g['m'],g['n'];G=gamma.reshape(m,n,-1);a=[];b=[];strength=[]
    for i in range(m+1):
        for j in range(n):
            a.append(X[i,j]);b.append(X[i,j+1]);strength.append((G[i,j] if i<m else 0)-(G[i-1,j] if i else 0))
    for i in range(m):
        for j in range(n+1):
            a.append(X[i,j]);b.append(X[i+1,j]);strength.append((G[i,j-1] if j else 0)-(G[i,j] if j<n else 0))
    K=kernel(P,np.array(a),np.array(b),eps);B=np.einsum('pij,pj->pi',K,N);s=np.array(strength);out=np.empty((len(P),gamma.shape[1]),complex)
    for row in range(len(P)):
        for col in range(gamma.shape[1]):out[row,col]=complex(math.fsum(B[row]*s[:,col].real),math.fsum(B[row]*s[:,col].imag))
    return out
def main():
    m=int(sys.argv[1]);assert m==128
    try:wrapper.experiment.run(m,192,1e-6)
    except AssertionError as err:
        tb=err.__traceback__;frame=None
        while tb:
            if 'sample' in tb.tb_frame.f_locals:frame=tb.tb_frame
            tb=tb.tb_next
        if frame is None:raise
        x=frame.f_locals;g=x['g'];sel=x['selected'];full=x['full'];U=unique_field(g,x['tail'],full,g['P'][sel],g['normal'][sel],x['eps']);R=U+x['WW']@full[-192:]-x['rr'];physical=np.linalg.norm(R,axis=0)/np.linalg.norm(x['rr'],axis=0);out=ROOT/'build/analysis/high-frequency-final/symmetry-diagnostic-m128.npz';np.savez_compressed(out,gamma=full,selected=sel,unique_residual=R,old_sample=x['sample']);result=dict(passed=bool(max(physical)<1e-11),scope='uncondensed physical unique-edge residual versus rounded ring-column diagnostic',ring_sample_residuals=x['sample'].tolist(),unique_edge_residuals=physical.tolist(),sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'rd/aero/unsteady-reference/experiments/run_high_frequency_final_bounded.py',ROOT/'rd/aero/unsteady-reference/experiments/run_high_frequency_final_symmetric.py']},artifacts={out.relative_to(ROOT).as_posix():hashlib.sha256(out.read_bytes()).hexdigest()});dest=ROOT/'rd/aero/unsteady-reference/records/high-frequency-symmetry-residual-diagnostic.json';assert not dest.exists();dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
    else:raise AssertionError('previous_failure_not_reproduced')
if __name__=='__main__':main()
