"""Synthetic policy branches and actual termination of a disposable child process."""
from pathlib import Path
import json,subprocess,sys,time
import avl_reference_resources as resource
import avl_span1280_watchdog as guard
ROOT=Path(__file__).resolve().parents[2]
def main():
    out=ROOT/'docs/verification/avl-span1280-watchdog-result.json'
    if out.exists():raise FileExistsError(out)
    good=dict(cgroup_known=True,effective_available=22*guard.GIB,swap_used=0)
    args=dict(elapsed=0,stage_elapsed=0,rss=0,console_bytes=0,total_bytes=0,snapshot=good,case_swap=0,stage_swap=0)
    cases=[]
    for name,changes,expected in [('normal',{},None),('case_time',dict(elapsed=1201),'per_process_time_limit'),
        ('stage_time',dict(stage_elapsed=2401),'stage_time_limit'),('rss',dict(rss=17*guard.GIB),'process_group_rss_limit'),
        ('console',dict(console_bytes=9*1024**2),'console_limit'),('disk',dict(total_bytes=513*1024**2),'disk_limit'),
        ('free',dict(snapshot=dict(good,effective_available=guard.GIB)),'available_floor'),
        ('swap',dict(snapshot=dict(good,swap_used=257*1024**2)),'swap_increase_limit'),
        ('cgroup',dict(snapshot=dict(good,cgroup_known=False)),'unknown_cgroup')]:
        actual=guard.stop_reason(**dict(args,**changes));assert actual==expected;cases.append(dict(name=name,reason=actual,passed=True))
    assert guard.start_reason(good) is None
    assert guard.start_reason(dict(good,cgroup_known=False))=='unknown_cgroup'
    assert guard.start_reason(dict(good,effective_available=19*guard.GIB))=='start_available_limit'
    proc=subprocess.Popen([sys.executable,'-c','import time; time.sleep(30)'],start_new_session=True)
    begin=time.monotonic()
    try:
        while time.monotonic()-begin<.3:time.sleep(.02)
        assert guard.stop_reason(**dict(args,elapsed=1201))=='per_process_time_limit';guard.terminate(proc)
        assert proc.returncode<0 and time.monotonic()-begin<5
    finally:guard.terminate(proc)
    result=dict(source_sha256=resource.digest(__file__),watchdog_sha256=resource.digest(guard.__file__),
        cases=cases,start_cases=3,actual_group_terminated=True,child_return_code=proc.returncode,seconds=time.monotonic()-begin,all_criteria_pass=True)
    out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
