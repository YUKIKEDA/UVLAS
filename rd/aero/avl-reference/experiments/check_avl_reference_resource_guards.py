"""Synthetic cgroup/memory observations: verify finite ancestor limits and unknown states."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
from unittest.mock import patch
import json
import avl_reference_resources as resource

ROOT=rd_paths.ROOT;GIB=1024**3

def main():
    records=[]
    specifications=[
        ('unlimited',None,1,None,2,True,22,True),
        ('finite_parent',None,1,20,5,True,15,False),
        ('finite_child',21,2,None,2,True,19,False),
        ('stricter_parent',23,1,20,5,True,15,False),
        ('exhausted_child',1,2,None,2,True,0,False),
        ('missing_controller',None,1,None,2,False,22,False)]
    for name,child_max,child_used,parent_max,parent_used,exposed,expected,allow in specifications:
        files={'/proc/meminfo':f'MemTotal: {24*GIB//1024} kB\nMemAvailable: {22*GIB//1024} kB\nSwapTotal: {6*GIB//1024} kB\nSwapFree: {6*GIB//1024} kB\n',
            '/proc/self/cgroup':'0::/job\n'}
        if exposed:
            for folder,maximum,used in [('/sys/fs/cgroup/job',child_max,child_used),('/sys/fs/cgroup',parent_max,parent_used)]:
                files[folder+'/memory.max']='max' if maximum is None else str(maximum*GIB)
                files[folder+'/memory.current']=str(used*GIB)
        def read(path,*args,**kwargs):
            key=path.as_posix()
            if key not in files:raise FileNotFoundError(key)
            return files[key]
        with patch.object(Path,'read_text',read),patch.object(Path,'exists',lambda p:p.as_posix() in files):
            observed=resource.memory_snapshot()
        eligible=observed['cgroup_known'] and observed['effective_available']>=20*GIB
        assert observed['effective_available']==expected*GIB and observed['cgroup_known']==exposed and eligible==allow
        records.append(dict(name=name,effective_available=observed['effective_available'],cgroup_known=observed['cgroup_known'],fine_start_eligible=eligible,passed=True))
    result=dict(source_sha256=resource.digest(__file__),resource_source_sha256=resource.digest(resource.__file__),
        cases=records,criteria_passed=len(records),all_criteria_pass=True,no_solver_run=True)
    (ROOT/'rd/aero/avl-reference/records/avl-span1280-resource-guard-result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
