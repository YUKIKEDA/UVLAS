"""Stage-specific resource policy for the read-only symmetric solve probe."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from avl_span1280_watchdog import terminate
GIB=1024**3
LIMITS={'small':(60,180,2*GIB),'main':(300,2400,5*GIB)}
def start_reason(snapshot):
    if not snapshot['cgroup_known']:return 'unknown_cgroup'
    if snapshot['effective_available']<8*GIB:return 'start_available_limit'
    return None
def stop_reason(stage,elapsed,total_elapsed,rss,console_bytes,total_bytes,snapshot,case_swap,stage_swap):
    per,total,mem=LIMITS[stage]
    if not snapshot['cgroup_known']:return 'unknown_cgroup'
    if elapsed>per:return 'per_process_time_limit'
    if total_elapsed>total:return 'stage_time_limit'
    if rss>mem:return 'process_group_rss_limit'
    if console_bytes>8*1024**2:return 'console_limit'
    if total_bytes>128*1024**2:return 'disk_limit'
    if snapshot['effective_available']<2*GIB:return 'available_floor'
    if snapshot['swap_used']-min(case_swap,stage_swap)>256*1024**2:return 'swap_increase_limit'
    return None
