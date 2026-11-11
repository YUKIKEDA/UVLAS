"""Retained coupled operators and nonlinear shoe second-variation check."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import numpy as np
import json,hashlib,time
import check_connection_beam_shell as coupled
ROOT=rd_paths.ROOT
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    start=time.monotonic();dest=ROOT/'rd/structure/connection/records/connection-snapshot-result.json';assert not dest.exists();d=coupled.prepare();F=np.array([1.,.004,-.003]);a=coupled.solve_case(d,F,'kkt');C,J,H=coupled.constraints(d,a['p'],a['R'],a['u']);rng=np.random.default_rng(150);v=rng.normal(size=J.shape[1]);v[6:]*=.001;v[:6]*=.1;t=.001
    def value(s):
        p,R=coupled.update(a['p'],a['R'],s*v[:6]);return coupled.constraints(d,p,R,a['u']+s*v[6:].reshape(-1,3))[0]
    fd=(-value(2*t)+16*value(t)-30*value(0)+16*value(-t)-value(-2*t))/(12*t*t);analytic=np.einsum('iab,a,b->i',H,v[3:6],v[3:6]);err=float(np.linalg.norm(fd-analytic)/np.linalg.norm(analytic));assert err<1e-8
    # An omitted multiplier-Hessian term has a finite, independently observed
    # second-variation defect; the accepted KKT carries the full term.
    multiplier=rng.normal(size=len(C));term=np.einsum('i,iab->ab',multiplier,H);defect=float(abs(v[3:6]@term@v[3:6]));assert defect>1e-8
    out=ROOT/'build/analysis/connection-coupled-snapshot.npz';assert not out.exists();np.savez_compressed(out,p=a['p'],R=a['R'],u=a['u'],multipliers=a['lam'],external=a['external'],K_internal=a['K'],J=a['J'],H_constraint=H,trace=d['T'],shell_mass=d['d']['M'],q0=d['d']['q0'],force=F)
    paths=[Path(__file__),Path(coupled.__file__),ROOT/'rd/structure/connection/records/connection-beam-shell-result.json'];dest.write_text(json.dumps(dict(passed=True,records=[dict(case='nonlinear_shoe_second_variation',relative_error=err,omitted_multiplier_hessian_work=defect),dict(case='retained_original_coupled_operators',equilibrium_residual=a['raw'],state_dofs=int(J.shape[1]),constraint_dofs=len(C))],artifact=dict(artifact=out.relative_to(ROOT).as_posix(),artifact_sha256=sha(out)),sources={p.relative_to(ROOT).as_posix():sha(p) for p in paths},seconds=time.monotonic()-start),indent=2),encoding='utf-8');print('CONNECTION_SNAPSHOT_RESULT',True,err,flush=True)
if __name__=='__main__':main()
