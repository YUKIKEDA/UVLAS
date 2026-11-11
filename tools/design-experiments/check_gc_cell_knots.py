"""Cheap rejection screen, before repeating the complete qualification."""
from pathlib import Path
import hashlib,json,time
import geometry_compatible_ancf_cell_knots as variant
import geometry_compatible_ancf_validation_stream as v
f=variant.f;ROOT=Path(__file__).resolve().parents[2]

def main():
    dest=ROOT/'docs/verification/geometry-compatible-ancf-cell-knots-screen.json';assert not dest.exists()
    records=[];geo=None
    for h in [.01,.0001]:
        tick=time.monotonic();v.select('rational');a=f.initialize(.4,4,0.,h,13,7,geo);geo=a['geo'];s=f.solve(a)
        o=f.observe(a,s,15,9) if s['passed_rank'] else None
        records.append(dict(h=h,observation=o,rank=s['passed_rank'],rcond=s.get('rcond'),seconds=time.monotonic()-tick,sizes=[a['nr'],a['nd'],a['nc']]))
        print(records[-1],flush=True)
    ratios={field:[y/max(x,1e-300) for x,y in zip(records[0]['observation'][field],records[1]['observation'][field])] for field in ['error','position']}
    paths=[Path(__file__),Path(variant.__file__)]
    result=dict(records=records,ratios=ratios,screen_passed=all(z<=2 for r in ratios.values() for z in r),main_shell_adopted=False,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(result['ratios'],result['screen_passed'],flush=True)

if __name__=='__main__':main()
