from pathlib import Path
import subprocess,time,json
root=Path.cwd();out=root/'out/issue108-contract';build=root/'out/issue94-build';results=[]
for config in ['Debug','Release']:
    steps=[('build',['cmake','--build',str(build),'--config',config,'--target','uvlas_test_structure_pg_coefficients','uvlas_test_structure_pg_coefficients_resources','uvlas_cli','-j','2'],600),('tests',['ctest','--test-dir',str(build),'-C',config,'--output-on-failure','-R',r'^structure\.pg_coefficients(_resources)?$'],650),('build-info',[str(build/config/'uvlas.exe'),'--build-info'],30)]
    for kind,command,limit in steps:
        logpath=out/('clean-'+config.lower()+'-'+kind+'.log');start=time.time()
        with logpath.open('wb') as log:
            r=subprocess.run(command,cwd=root,stdout=log,stderr=subprocess.STDOUT,timeout=limit)
        row={'configuration':config,'step':kind,'command':command,'exit_code':r.returncode,'seconds':time.time()-start,'log':str(logpath)};results.append(row)
        (out/'clean-validation-results.json').write_text(json.dumps(results,indent=2))
        print(config,kind,r.returncode,round(row['seconds'],3),flush=True)
        if r.returncode:raise SystemExit(r.returncode)
