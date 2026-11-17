from pathlib import Path
import json,subprocess,time
root=Path.cwd();out=root/'out/issue108-tidy';results=[]
for name in ['src/structure/pg_coefficients.cpp','tests/structure/pg_coefficients_test.cpp']:
    path=root/name;start=time.time();logpath=out/(path.stem+'-final.log')
    with logpath.open('wb') as log:
        r=subprocess.run(['uvx','--offline','clang-tidy@21.1.1','-p',str(out),'--quiet',str(path)],cwd=root,stdout=log,stderr=subprocess.STDOUT,timeout=240)
    results.append({'file':name,'exit_code':r.returncode,'seconds':time.time()-start,'log':str(logpath)})
    print(name,r.returncode,flush=True)
(out/'final-results.json').write_text(json.dumps(results,indent=2))
