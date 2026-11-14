import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import sys,json,hashlib,numpy as np
from scipy.linalg import solve
from uvlm_frequency_delta import Solver
from frequency_motion_columns import observe
ROOT=rd_paths.ROOT
def main():
    side=sys.argv[1];records=[]
    for m,n in [(2,8),(4,24)]:
        if side=='candidate':
            from uvlm_frequency_fast import geometry,bound,age_nodes,age_integral
            g=geometry(m,n)
        else:
            from ptera_frequency_reference import geometry,run,project
            g=geometry(m,n)
        s=Solver(g,side)
        for k in [.05,.2,.8]:
            w=40*k/1.114285714285714;r=s.run(k);X=g['X']
            if side=='candidate':
                A=bound(g['P'],g['normal'],X,1e-5).astype(complex);ages,weights=age_nodes(g,w,128.,16);W=age_integral(g['P'],g['normal'],X[-1],1e-5,ages,weights)
            else:
                old=run(m,n,k,1e-5,128.,16,True);W=old['W'];A=project(g,X[:-1,:-1].reshape(-1,3),X[:-1,1:].reshape(-1,3),X[1:,1:].reshape(-1,3),X[1:,:-1].reshape(-1,3),1e-5).astype(complex)
            A[:,-n:]+=W;rhs=np.column_stack([1j*w*g['normal'][:,2],-(20+1j*w*(g['P'][:,0]-.375))*g['normal'][:,2]]);gamma=solve(A,rhs);transfer=observe(g,gamma,w,side);errors=dict(matrix=float(np.linalg.norm(A-(s.A.astype(complex)+np.pad(r['W'],((0,0),(m*n-n,0)))))/np.linalg.norm(A)),gamma=float(np.linalg.norm(gamma-r['gamma'])/np.linalg.norm(gamma)),transfer=float(np.max(abs((transfer-r['transfer'])/transfer))));assert max(errors.values())<1e-9;records.append(dict(m=m,n=n,k=k,errors=errors,residuals=r['residuals']))
    sources={}
    for p in [Path(__file__),ROOT/'rd/aero/unsteady-reference/experiments/uvlm_frequency_delta.py',ROOT/'rd/aero/unsteady-reference/records/frequency-delta-basis-protocol.md']:
        sources[p.relative_to(ROOT).as_posix()]=hashlib.sha256(p.read_bytes()).hexdigest()
    dest=ROOT/rd_paths.record(f'frequency-delta-{side}-result.json');assert not dest.exists();dest.write_text(json.dumps(dict(passed=True,records=records,sources=sources),indent=2),encoding='utf-8');print(json.dumps(dict(passed=True,side=side,max_error=max(max(r['errors'].values()) for r in records))),flush=True)
if __name__=='__main__':main()
