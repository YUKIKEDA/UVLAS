"""Reference-only entrance gates for the independently enriched trial spaces."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,importlib,json
ROOT=rd_paths.ROOT


def main():
    dest=ROOT/'rd/structure/main-shell/paired-geometry/records/pg-ritz-split-screen-result.json';assert not dest.exists()
    records=[];paths=[Path(__file__)]
    for name in ['pg_ritz_shell_v4','pg_ritz_shell_v5','pg_ritz_shell_v6',
                 'pg_ritz_shell_v7','pg_ritz_shell_v8','pg_ritz_shell_v9',
                 'pg_ritz_shell_v10','pg_ritz_shell_particular_only']:
        module=importlib.import_module(name);paths.extend([Path(module.__file__),Path(module.lift.__file__)])
        try:
            m=module.prepare(2,(.4,.2,0.));r=dict(module=name,entrance_passed=True,ns=m['ns'])
        except ValueError as exc:r=dict(module=name,entrance_passed=False,reason=str(exc))
        records.append(r);print(r,flush=True)
    assert all(not r['entrance_passed'] for r in records)
    result=dict(accepted=False,records=records,R4_complete=False,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    dest.write_text(json.dumps(result,indent=2),encoding='utf-8')


if __name__=='__main__':main()
