from pathlib import Path
import subprocess,time,json
root=Path.cwd();out=root/'out/issue108-evidence';build=root/'out/issue94-build';rows=[]
for config in ['Debug','Release']:
    for kind,command,limit in [('build',['cmake','--build',str(build),'--config',config,'--target','uvlas_test_structure_pg_coefficients','uvlas_cli','-j','2'],180),('jet',[str(build/'tests'/config/'uvlas_test_structure_pg_coefficients.exe'),'PG all quartic*'],120),('build-info',[str(build/config/'uvlas.exe'),'--build-info'],30)]:
        logpath=out/('jet-'+config.lower()+'-'+kind+'.log');start=time.time()
        with logpath.open('wb') as log:r=subprocess.run(command,cwd=root,stdout=log,stderr=subprocess.STDOUT,timeout=limit)
        rows.append({'configuration':config,'kind':kind,'command':command,'exit_code':r.returncode,'seconds':time.time()-start,'log':str(logpath)})
        (out/'jet-validation-results.json').write_text(json.dumps(rows,indent=2));print(config,kind,r.returncode,flush=True)
        if r.returncode:raise SystemExit(r.returncode)
with (out/'jet-normal-tidy.log').open('wb') as log:r=subprocess.run(['uvx','--offline','clang-tidy@21.1.1','-p',str(root/'out/issue108-tidy'),'--quiet',str(root/'tests/structure/pg_coefficients_test.cpp')],cwd=root,stdout=log,stderr=subprocess.STDOUT,timeout=120)
print('normal tidy',r.returncode,flush=True)
(out/'jet-tidy-result.json').write_text(json.dumps({'exit_code':r.returncode,'command':'uvx --offline clang-tidy@21.1.1 -p out/issue108-tidy --quiet tests/structure/pg_coefficients_test.cpp'}))
raise SystemExit(r.returncode)
