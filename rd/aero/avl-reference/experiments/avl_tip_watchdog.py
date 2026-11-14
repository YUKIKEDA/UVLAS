"""Fixed-count redistribution resource policy, separate from the N1280 fine policy."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from avl_span1280_watchdog import terminate
GIB=1024**3
def start_reason(snapshot):
    if not snapshot['cgroup_known']:return 'unknown_cgroup'
    if snapshot['effective_available']<8*GIB:return 'start_available_limit'
    return None
def stop_reason(elapsed,stage_elapsed,rss,console_bytes,total_bytes,snapshot,case_swap,stage_swap):
    if not snapshot['cgroup_known']:return 'unknown_cgroup'
    if elapsed>300:return 'per_process_time_limit'
    if stage_elapsed>1200:return 'stage_time_limit'
    if rss>5*GIB:return 'process_group_rss_limit'
    if console_bytes>8*1024**2:return 'console_limit'
    if total_bytes>128*1024**2:return 'disk_limit'
    if snapshot['effective_available']<2*GIB:return 'available_floor'
    if snapshot['swap_used']-min(case_swap,stage_swap)>256*1024**2:return 'swap_increase_limit'
    return None
