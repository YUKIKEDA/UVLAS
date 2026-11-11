"""Predeclared fine-run stop policy and process-group termination."""
import os,signal
GIB=1024**3
def start_reason(snapshot):
    if not snapshot['cgroup_known']:return 'unknown_cgroup'
    if snapshot['effective_available']<20*GIB:return 'start_available_limit'
    return None
def stop_reason(elapsed,stage_elapsed,rss,console_bytes,total_bytes,snapshot,case_swap,stage_swap):
    if not snapshot['cgroup_known']:return 'unknown_cgroup'
    if elapsed>1200:return 'per_process_time_limit'
    if stage_elapsed>2400:return 'stage_time_limit'
    if rss>16*GIB:return 'process_group_rss_limit'
    if console_bytes>8*1024**2:return 'console_limit'
    if total_bytes>512*1024**2:return 'disk_limit'
    if snapshot['effective_available']<2*GIB:return 'available_floor'
    if snapshot['swap_used']-min(case_swap,stage_swap)>256*1024**2:return 'swap_increase_limit'
    return None
def terminate(proc):
    if proc.poll() is None:
        try:os.killpg(proc.pid,signal.SIGKILL)
        except ProcessLookupError:pass
    proc.wait()
