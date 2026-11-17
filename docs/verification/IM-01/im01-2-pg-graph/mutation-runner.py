from pathlib import Path
import json,subprocess,hashlib,re,shutil
root=Path.cwd();out=root/'out/issue94-mutations';build=root/'out/issue94-build';snap=out/'snapshot';snap.mkdir(parents=True,exist_ok=True)
sha=lambda b:hashlib.sha256(b).hexdigest()
inputs={}
for rel in ['src/structure/pg_graph.cpp','src/structure/pg_graph_detail.hpp','src/structure/exact_polynomial_detail.hpp','tests/structure/pg_graph_test.cpp','tests/structure/fixtures/pg-graph-v1.json']:
 src=root/rel;dst=snap/rel;dst.parent.mkdir(parents=True,exist_ok=True);data=src.read_bytes();dst.write_bytes(data);assert dst.read_bytes()==data;inputs[rel]=sha(data)
for src in (root/'include').rglob('*'):
 if src.is_file():
  rel=src.relative_to(root);dst=snap/rel;dst.parent.mkdir(parents=True,exist_ok=True);data=src.read_bytes();dst.write_bytes(data);inputs[rel.as_posix()]=sha(data)
source=(snap/'src/structure/pg_graph.cpp').read_bytes()
db=json.loads((build/'compile_commands.json').read_text());entry=next(x for x in db if x['file'].endswith('pg_graph.cpp') and 'Release' in x['command'] and 'uvlas_structure.dir' in x['command'])
cmd=entry['command'];cmd=cmd.replace(str(root/'include').replace('\\','/'),str(snap/'include').replace('\\','/')).replace(str(root/'include'),str(snap/'include'))
ninja=(build/'CMakeFiles/impl-Release.ninja').read_text();part=ninja[ninja.index('build tests\\Release\\uvlas_test_structure_pg_graph.exe:'):];part=part[:part.index('\n\n')];libs=re.search(r'  LINK_LIBRARIES = (.*)',part)[1].split()
libs=[str(build/x) if x.startswith('Release\\') else x for x in libs]
testobj=build/'tests/CMakeFiles/uvlas_test_structure_pg_graph.dir/Release/structure/pg_graph_test.cpp.obj'
for p in [testobj,build/'Release/uvlas_structure.lib',build/'Release/uvlas_core.lib',build/'compile_commands.json']:
 inputs[p.relative_to(root).as_posix()]=sha(p.read_bytes())
cases=[('control',None,None),('lower_side',b'if (side == PgSide::Lower && boundary > 0)\n                return boundary - 1;',b'if (side == PgSide::Lower && boundary > 0)\n                return boundary;'),('local_scale',b'denominator = mul(denominator, static_cast<std::uint64_t>(4 * d.n));\n                const auto scale',b'denominator = mul(denominator, static_cast<std::uint64_t>(2 * d.n));\n                const auto scale'),('identity_skew',b'    append_rational(hash, d.skew);',b'    // mutant: omitted skew identity contribution'),('cell_index',b'd.cells.values()[static_cast<std::size_t>(a * axis + b)] = {a, b};',b'd.cells.values()[static_cast<std::size_t>(a * axis + b)] = {a, 0};')]
results=[]
for name,old,new in cases:
 if name in ['control','lower_side']:
  case=out/name;data=(case/'pg_graph.cpp').read_bytes()
  assert (case/'producer.obj').exists() and (case/'run.log').exists()
  results.append({'case':name,'compile_exit':0,'link_exit':0,'run_exit':0 if name=='control' else 42,'source_sha256':sha(data),'original_archive_member_extracted':False,'completed_by_prior_invocation':True,'old_bytes':None if old is None else old.decode(),'new_bytes':None if new is None else new.decode()})
  continue
 case=out/name;case.mkdir(exist_ok=True);src=case/'pg_graph.cpp';data=source
 if old is not None:assert data.count(old)==1,(name,data.count(old));data=data.replace(old,new)
 src.write_bytes(data)
 if name=='control':assert src.read_bytes()==(root/'src/structure/pg_graph.cpp').read_bytes()
 obj=case/'producer.obj';exe=case/'test.exe';mapping=case/'link.map'
 command=cmd.replace(entry['file'],str(src)).replace(entry['file'].replace('\\','/'),str(src))
 command=re.sub(r'/Fo\S+',lambda m:'/Fo'+str(obj),command)
 command+=' /I'+str(snap/'src/structure')
 (case/'compile-command.txt').write_text(command,encoding='utf-8')
 with (case/'compile.log').open('wb') as f:r=subprocess.run(command,cwd=entry['directory'],stdout=f,stderr=subprocess.STDOUT,timeout=120)
 assert r.returncode==0,(name,'compile',r.returncode)
 link=['link.exe','/NOLOGO','/INCREMENTAL:NO','/SUBSYSTEM:CONSOLE','/OUT:'+str(exe),'/MAP:'+str(mapping),str(obj),str(testobj),*libs]
 (case/'link-command.json').write_text(json.dumps(link),encoding='utf-8')
 with (case/'link.log').open('wb') as f:r=subprocess.run(link,cwd=build,stdout=f,stderr=subprocess.STDOUT,timeout=30)
 assert r.returncode==0,(name,'link',r.returncode)
 with (case/'run.log').open('wb') as f:r=subprocess.run([str(exe)],cwd=root,stdout=f,stderr=subprocess.STDOUT,timeout=30)
 maptext=mapping.read_text();assert 'uvlas_structure:pg_graph.cpp.obj' not in maptext
 assert r.returncode==(0 if name=='control' else 42),(name,r.returncode)
 results.append({'case':name,'compile_exit':0,'link_exit':0,'run_exit':r.returncode,'source_sha256':sha(data),'original_archive_member_extracted':False,'old_bytes':None if old is None else old.decode(),'new_bytes':None if new is None else new.decode()});print(name,r.returncode,flush=True)
for rel,digest in inputs.items():assert sha((root/rel).read_bytes())==digest,rel
(out/'results.json').write_text(json.dumps({'source_commit':'1ae6d19f91a95b820ce3f2bf17f8907cdf2d5f16','inputs':inputs,'results':results,'original_input_mismatches':0},indent=2),encoding='utf-8')
