"""Linux reference experiment resource observations; no solver configuration changes."""
from pathlib import Path
import hashlib, json, re, time

GIB=1024**3

def memory_snapshot():
    info={k:int(v)*1024 for k,v in re.findall(r'^(\w+):\s*(\d+) kB',Path('/proc/meminfo').read_text(),re.M)}
    controls=[];known=False;reason=None
    try:
        line=next(s for s in Path('/proc/self/cgroup').read_text().splitlines() if s.startswith('0::'))
        relative=line.split('::',1)[1].lstrip('/')
        # Require the exposed cgroup path, rather than guessing a different namespace root.
        root=Path('/sys/fs/cgroup');folder=root/relative
        assert folder.is_relative_to(root) and (folder/'memory.max').exists()
        while True:
            maximum=(folder/'memory.max').read_text().strip();current=int((folder/'memory.current').read_text())
            controls.append(dict(path=str(folder),maximum=None if maximum=='max' else int(maximum),current=current))
            if folder==root:break
            folder=folder.parent
        known=True
    except (OSError,ValueError,StopIteration,AssertionError) as exc:reason=repr(exc)
    limits=[max(0,r['maximum']-r['current']) for r in controls if r['maximum'] is not None]
    return dict(mem_total=info['MemTotal'],mem_available=info['MemAvailable'],
        swap_used=info.get('SwapTotal',0)-info.get('SwapFree',0),cgroup_known=known,cgroup_unknown_reason=reason,
        cgroups=controls,effective_available=min([info['MemAvailable'],*limits]))

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
    import ctypes.util,os,platform,shutil
    root=Path(__file__).resolve().parents[2];folder=root/'docs/verification'
    output=folder/'avl-span1280-resource-audit.json'
    if output.exists():raise FileExistsError(output)
    base=root/'build/reference-avl/source-symmetry8000-spacing-pi64/src'
    previous=json.loads((folder/'avl-symmetry-pi-fine-runs.json').read_text())
    large=[r for r in previous['runs'] if r['m']==24 and r['physical_n']==640]
    assert len(large)==2
    matrix=16000**2*8*7
    result=dict(source_sha256=digest(__file__),protocol_sha256=digest(folder/'avl-span1280-capacity-experiment.md'),
        observed_unix=time.time(),platform=platform.platform(),logical_cpus=os.cpu_count(),memory=memory_snapshot(),
        disk_free=shutil.disk_usage(root).free,external_libraries={n:ctypes.util.find_library(n) for n in ['blas','lapack','openblas']},
        source_files={n:digest(base/n) for n in ['avl_heap.f90','avl_heap_inc.f90','ADIMEN.INC','matrix-lapackdp.f','matrix-lapacksubs.f','asetup.f']},
        parent_runs_sha256=digest(folder/'avl-symmetry-pi-fine-runs.json'),parent_max_seconds=max(r['seconds'] for r in large),
        parent_peak_bytes=max(r['peak_resident_bytes'] for r in large),allocated_dense_bytes=matrix,allocated_dense_gib=matrix/GIB,
        active_dense_bytes=15360**2*8*7,cubic_time_scenario=max(r['seconds'] for r in large)*8,
        allocation_is_not_peak_bound=True,time_is_not_performance_guarantee=True,
        proposed_limits=dict(regression_rss=5*GIB,fine_rss=16*GIB,fine_start_available=20*GIB,fine_remaining_floor=2*GIB),
        no_fine_run_performed=True)
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['memory','allocated_dense_bytes','allocated_dense_gib','cubic_time_scenario','external_libraries']},indent=2))

if __name__=='__main__':main()
