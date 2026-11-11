"""Actual >100k coordinate storage/action/solve using tested two-patch blocks.

This is a bounded-patch capacity experiment, not a claim about a 100k globally
coupled production solve or product runtime acceptance.
"""
from pathlib import Path
import ctypes,hashlib,json,time
import numpy as np
from scipy.linalg import solve_triangular
ROOT=Path(__file__).resolve().parents[2]

def memory():
    class Counts(ctypes.Structure):
        _fields_=[('cb',ctypes.c_ulong),('PageFaultCount',ctypes.c_ulong)]+[(k,ctypes.c_size_t) for k in
            ['PeakWorkingSetSize','WorkingSetSize','QuotaPeakPagedPoolUsage','QuotaPagedPoolUsage','QuotaPeakNonPagedPoolUsage','QuotaNonPagedPoolUsage','PagefileUsage','PeakPagefileUsage']]
    out=Counts();out.cb=ctypes.sizeof(out)
    kernel=ctypes.WinDLL('kernel32');kernel.GetCurrentProcess.restype=ctypes.c_void_p
    psapi=ctypes.WinDLL('psapi');psapi.GetProcessMemoryInfo.argtypes=[ctypes.c_void_p,ctypes.POINTER(Counts),ctypes.c_ulong]
    if not psapi.GetProcessMemoryInfo(kernel.GetCurrentProcess(),ctypes.byref(out),out.cb):raise OSError('GetProcessMemoryInfo')
    return dict(rss=out.WorkingSetSize,peak_rss=out.PeakWorkingSetSize)

def main():
    dest=ROOT/'docs/verification/mass-paired-capacity-result.json';assert not dest.exists();start=time.monotonic()
    parent=ROOT/'docs/verification/mass-paired-interface-result.json';report=json.loads(parent.read_text());assert report['passed']
    path=ROOT/report['records'][1]['artifact'];assert hashlib.sha256(path.read_bytes()).hexdigest()==report['records'][1]['artifact_sha256']
    a=np.load(path);T=a['root'];scale=a['scale'];dimension=len(scale);count=(100000+dimension-1)//dimension
    before=memory();tick=time.monotonic()
    # Each copy occupies actual memory: no broadcast-only capacity estimate.
    roots=np.repeat(T[None],count,axis=0);scales=np.repeat(scale[None],count,axis=0)
    rng=np.random.default_rng(20261105);x=rng.normal(size=(count,dimension))*scales
    allocation=time.monotonic()-tick;tick=time.monotonic()
    y=np.einsum('bij,bj->bi',roots,x/scales,optimize=True)
    rhs=np.einsum('bji,bj->bi',roots,y,optimize=True)/scales
    action=time.monotonic()-tick;tick=time.monotonic();solved=np.empty_like(x)
    for b in range(count):
        first=solve_triangular(roots[b].T,scales[b]*rhs[b],lower=True,check_finite=False)
        solved[b]=scales[b]*solve_triangular(roots[b],first,check_finite=False)
    solve=time.monotonic()-tick
    error=float(np.linalg.norm(np.einsum('bij,bj->bi',roots,(solved-x)/scales))/np.linalg.norm(y))
    residual=np.einsum('bij,bj->bi',roots,solved/scales,optimize=True)
    rr=np.einsum('bji,bj->bi',roots,residual,optimize=True)-scales*rhs
    backward=float(np.linalg.norm(rr)/np.linalg.norm(scales*rhs));after=memory()
    checks=dict(actual_dofs=count*dimension>=100000,energy_solution=error<=1e-8,backward=backward<=1e-8,
        memory=after['peak_rss']<4*1024**3,time=time.monotonic()-start<=300.)
    # Constant local dense factors and explicit face restriction: no global dense N/M/Z.
    record=dict(blocks=count,block_dofs=dimension,total_dofs=count*dimension,root_bytes=roots.nbytes,
        state_and_scale_bytes=x.nbytes+solved.nbytes+scales.nbytes,allocation_seconds=allocation,action_seconds=action,solve_seconds=solve,
        energy_error=error,backward_error=backward,memory_before=before,memory_after=after,checks=checks,
        scope='independent copies of the verified nonmatching two-curved-patch supported physical problem; connected 100k solver performance not measured')
    result=dict(passed=all(checks.values()),record=record,seconds=time.monotonic()-start,R4_complete=False,
        sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),parent,path]})
    dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
