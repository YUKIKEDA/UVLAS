from pathlib import Path
import json,re,subprocess,hashlib
root=Path.cwd();build=root/'out/issue94-build';out=root/'out/issue108-contract';source=out/'public-n4-measure.cpp';obj=out/'public-n4-measure.obj';exe=out/'public-n4-measure.exe'
entry=next(x for x in json.loads((build/'compile_commands.json').read_text()) if x['file'].endswith('pg_coefficients_test.cpp') and 'Release' in x['command'])
command=entry['command'].replace(entry['file'],str(source)).replace(entry['file'].replace('\\','/'),str(source))
command=re.sub(r'/Fo\S+',lambda m:'/Fo'+str(obj),command)
(out/'measure-compile-command.txt').write_text(command)
with (out/'measure-compile.log').open('wb') as log:r=subprocess.run(command,cwd=entry['directory'],stdout=log,stderr=subprocess.STDOUT,timeout=120)
assert r.returncode==0,('compile',r.returncode)
ninja=(build/'CMakeFiles/impl-Release.ninja').read_text();part=ninja[ninja.index('build tests\\Release\\uvlas_test_structure_pg_coefficients.exe:'):];part=part[:part.index('\n\n')];libs=re.search(r'  LINK_LIBRARIES = (.*)',part)[1].split();libs=[str(build/x) if x.startswith('Release\\') else x for x in libs]
command=['link.exe','/NOLOGO','/INCREMENTAL:NO','/SUBSYSTEM:CONSOLE','/OUT:'+str(exe),str(obj),*libs]
(out/'measure-link-command.json').write_text(json.dumps(command))
with (out/'measure-link.log').open('wb') as log:r=subprocess.run(command,cwd=build,stdout=log,stderr=subprocess.STDOUT,timeout=30)
assert r.returncode==0,('link',r.returncode)
print('measurement harness compile/link0 sha256='+hashlib.sha256(source.read_bytes()).hexdigest())
