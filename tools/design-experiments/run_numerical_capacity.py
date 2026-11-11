"""Bounded capacity witness; neither synthetic matrix is a physical PG/UVLM case."""
from pathlib import Path
import subprocess,sys,time,json,platform,hashlib
ROOT=Path(__file__).resolve().parents[2]
def main():
    executable=Path(sys.argv[1]).resolve();name=sys.argv[2];out=ROOT/f'build/analysis/numerical-capacity-{name}.json';assert not out.exists()
    before=time.perf_counter()
    p=subprocess.run([str(executable)],capture_output=True,text=True,timeout=600,check=False)
    log=ROOT/f'build/analysis/numerical-capacity-{name}.log';log.write_text(p.stdout+'\nSTDERR\n'+p.stderr,encoding='utf-8')
    values=[json.loads(line) for line in p.stdout.splitlines() if line.startswith('{')]
    result=dict(passed=p.returncode==0 and len(values)==2,scope='capacity and algebra only; not physical assembly or product performance acceptance',platform=platform.platform(),elapsed_seconds=time.perf_counter()-before,executable_sha256=hashlib.sha256(executable.read_bytes()).hexdigest(),cases=values)
    out.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result));assert result['passed']
if __name__=='__main__':main()
