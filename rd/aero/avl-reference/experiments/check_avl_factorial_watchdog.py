"""Policy normal/boundary/refusal cases and real process-group termination."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import json,subprocess,sys,time
import avl_reference_resources as resource
import avl_factorial_watchdog as guard
import avl_span1280_watchdog as termination
ROOT=rd_paths.ROOT
def main():
    output=ROOT/'rd/aero/avl-reference/records/avl-factorial-watchdog-result.json'
    if output.exists():raise FileExistsError(output)
    good=dict(cgroup_known=True,effective_available=10*guard.GIB,swap_used=0)
    args=dict(elapsed=0,stage_elapsed=0,rss=0,console_bytes=0,total_bytes=0,snapshot=good,case_swap=0,stage_swap=0)
    cases=[]
    specs=[('normal',{},None),('boundary',dict(elapsed=400,stage_elapsed=1600,rss=6*guard.GIB,console_bytes=8*1024**2,total_bytes=128*1024**2,snapshot=dict(good,effective_available=2*guard.GIB,swap_used=256*1024**2)),None),
        ('time',dict(elapsed=401),'per_process_time_limit'),('stage',dict(stage_elapsed=1601),'stage_time_limit'),('rss',dict(rss=6*guard.GIB+1),'process_group_rss_limit'),
        ('console',dict(console_bytes=8*1024**2+1),'console_limit'),('disk',dict(total_bytes=128*1024**2+1),'disk_limit'),
        ('available',dict(snapshot=dict(good,effective_available=2*guard.GIB-1)),'available_floor'),('swap',dict(snapshot=dict(good,swap_used=256*1024**2+1)),'swap_increase_limit'),
        ('cgroup',dict(snapshot=dict(good,cgroup_known=False)),'unknown_cgroup')]
    for name,changes,expected in specs:
        actual=guard.stop_reason(**dict(args,**changes));assert actual==expected;cases.append(dict(name=name,reason=actual,passed=True))
    assert guard.start_reason(good) is None
    assert guard.start_reason(dict(good,effective_available=10*guard.GIB-1))=='start_available_limit'
    assert guard.start_reason(dict(good,cgroup_known=False))=='unknown_cgroup'
    proc=subprocess.Popen([sys.executable,'-c','import time;time.sleep(30)'],start_new_session=True);begin=time.monotonic()
    try:
        time.sleep(.3);assert guard.stop_reason(**dict(args,elapsed=401))=='per_process_time_limit';guard.terminate(proc);assert proc.returncode<0 and time.monotonic()-begin<5
    finally:guard.terminate(proc)
    result=dict(source_sha256=resource.digest(__file__),watchdog_sha256=resource.digest(guard.__file__),termination_source_sha256=resource.digest(termination.__file__),
        cases=cases,start_cases=3,real_termination=True,child_return_code=proc.returncode,all_criteria_pass=True)
    output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(cases=10,start_cases=3,real_termination=True,all_criteria_pass=True)))
if __name__=='__main__':main()
