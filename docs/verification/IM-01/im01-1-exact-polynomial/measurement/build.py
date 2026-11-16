from pathlib import Path
import json,re,subprocess
root=Path.cwd();out=root/'out/issue87-measure';build=root/'out/issue87-build';db=json.loads((build/'compile_commands.json').read_text());entry=next(x for x in db if x['file'].endswith('exact_polynomial.cpp') and 'Release' in x['command'] and 'uvlas_structure.dir' in x['command']);cmd=entry['command'].replace(entry['file'],str(out/'probe.cpp'));cmd=re.sub(r'/Fo\S+',lambda m:'/Fo'+str(out/'probe.obj'),cmd);cmd+=' /I'+str(root/'src/structure');(out/'compile-command.txt').write_text(cmd)
with (out/'compile.log').open('wb') as f:r=subprocess.run(cmd,cwd=entry['directory'],stdout=f,stderr=subprocess.STDOUT,timeout=120)
assert r.returncode==0
link=['link.exe','/NOLOGO','/INCREMENTAL:NO','/OUT:'+str(out/'probe.exe'),str(out/'probe.obj'),str(build/'Release/uvlas_structure.lib'),str(build/'Release/uvlas_core.lib'),'psapi.lib'];(out/'link-command.json').write_text(json.dumps(link))
with (out/'link.log').open('wb') as f:r=subprocess.run(link,cwd=root,stdout=f,stderr=subprocess.STDOUT,timeout=30)
assert r.returncode==0
