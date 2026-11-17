from pathlib import Path
import subprocess,json,re
root=Path.cwd();out=root/'out/issue87-evidence';build=root/'out/issue87-build'
assert not subprocess.check_output(['git','status','--porcelain','--untracked-files=no'],cwd=root).strip()
def run(name,cmd,timeout=180):
 with (out/name).open('wb') as f:r=subprocess.run(cmd,cwd=root,stdout=f,stderr=subprocess.STDOUT,timeout=timeout)
 assert r.returncode==0,(name,r.returncode);print(name,r.returncode,flush=True)
db=json.loads((build/'compile_commands.json').read_text());entry=next(x for x in db if x['file'].endswith('exact_polynomial.cpp') and 'Debug' in x['command'] and 'uvlas_structure.dir' in x['command']);src=out/'build-info.cpp';src.write_text('#include <uvlas/core/build_info.hpp>\n#include <iostream>\nint main(){std::cout<<uvlas::build_info_json()<<"\\n";}\n');cmd=entry['command'].replace(entry['file'],str(src));cmd=re.sub(r'/Fo\S+',lambda m:'/Fo'+str(out/'build-info.obj'),cmd);cmd=re.sub(r'/Fd\S+',lambda m:'/Fd'+str(out/'build-info.pdb'),cmd);run('build-info-compile.log',cmd)
run('build-info-link.log',['link.exe','/NOLOGO','/OUT:'+str(out/'build-info.exe'),str(out/'build-info.obj'),str(build/'Debug/uvlas_core.lib'),'psapi.lib'])
run('debug-build-info.json',[str(out/'build-info.exe')]);info=json.loads((out/'debug-build-info.json').read_text());assert info['uncommitted_changes'] is False
print('clean debug source',info.get('commit'),flush=True)
run('measure-build.log',['uv','run','--no-project','python',str(root/'out/issue87-measure/build.py')])
