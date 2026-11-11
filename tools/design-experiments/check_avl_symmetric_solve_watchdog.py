"""Check normal/boundary/refusal branches for both probe stages and real termination."""
from pathlib import Path
import json,subprocess,sys,time
import avl_symmetric_solve_watchdog as guard
import avl_span1280_watchdog as termination
import avl_reference_resources as resource
ROOT=Path(__file__).resolve().parents[2]
def main():
    output=ROOT/'docs/verification/avl-symmetric-solve-watchdog-result.json'
    if output.exists():raise FileExistsError(output)
    good=dict(cgroup_known=True,effective_available=8*guard.GIB,swap_used=0);cases=[]
    for stage,(per,total,rss) in guard.LIMITS.items():
        base=dict(stage=stage,elapsed=0,total_elapsed=0,rss=0,console_bytes=0,total_bytes=0,snapshot=good,case_swap=0,stage_swap=0)
        specs=[('normal',{},None),('boundary',dict(elapsed=per,total_elapsed=total,rss=rss,console_bytes=8*1024**2,total_bytes=128*1024**2,snapshot=dict(good,effective_available=2*guard.GIB,swap_used=256*1024**2)),None),
            ('time',dict(elapsed=per+1),'per_process_time_limit'),('stage',dict(total_elapsed=total+1),'stage_time_limit'),('rss',dict(rss=rss+1),'process_group_rss_limit'),
            ('console',dict(console_bytes=8*1024**2+1),'console_limit'),('disk',dict(total_bytes=128*1024**2+1),'disk_limit'),
            ('floor',dict(snapshot=dict(good,effective_available=2*guard.GIB-1)),'available_floor'),('swap',dict(snapshot=dict(good,swap_used=256*1024**2+1)),'swap_increase_limit'),
            ('cgroup',dict(snapshot=dict(good,cgroup_known=False)),'unknown_cgroup')]
        for label,changes,expected in specs:
            assert guard.stop_reason(**dict(base,**changes))==expected;cases.append(dict(stage=stage,case=label,reason=expected,passed=True))
    assert guard.start_reason(good) is None and guard.start_reason(dict(good,effective_available=8*guard.GIB-1))=='start_available_limit' and guard.start_reason(dict(good,cgroup_known=False))=='unknown_cgroup'
    proc=subprocess.Popen([sys.executable,'-c','import time;time.sleep(30)'],start_new_session=True);begin=time.monotonic()
    try:time.sleep(.3);guard.terminate(proc);assert proc.returncode<0 and time.monotonic()-begin<5
    finally:guard.terminate(proc)
    result=dict(source_sha256=resource.digest(__file__),watchdog_sha256=resource.digest(guard.__file__),termination_sha256=resource.digest(termination.__file__),cases=cases,start_cases=3,real_termination=True,all_criteria_pass=True)
    output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(cases=len(cases),real_termination=True,all_criteria_pass=True)))
if __name__=='__main__':main()
