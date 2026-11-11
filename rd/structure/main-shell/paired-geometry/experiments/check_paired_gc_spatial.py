"""Complete non-floor spatial evidence, using the already fixed basis and targets."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import gc,hashlib,json,time
import check_paired_gc_qualification as full
v=full.v;f=v.f;ROOT=full.ROOT

def main():
    parent=ROOT/'rd/structure/main-shell/paired-geometry/records/paired-geometry-ancf-candidate-result.json';base=json.loads(parent.read_text());assert base['adoption_candidate_qualified']
    dest=ROOT/'rd/structure/main-shell/paired-geometry/records/paired-geometry-ancf-spatial-result.json';assert not dest.exists();start=time.monotonic();records=[];stops=[];geos={};out=ROOT/'build/analysis'/('paired-geometry-ancf-spatial-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir()
    for amp in [0.,.2]:
        for h in [.01,.0001]:
            try:
                tick=time.monotonic();v.select('rational');a=f.initialize(.4,2,amp,h,13,7,geos.get(amp));geos[amp]=a['geo'];s=f.solve(a);assert s['passed_rank'];primary=f.observe(a,s,13,7);ind=f.observe(a,s,15,9);traces=v.interface(a,s)
                checks=dict(rank=s['passed_rank'],condition=s['rcond']>=1e-10,linear=max(s['linear'])<=1e-8,original_stationarity=max(primary['original_kkt'])<=1e-8,independent=max(abs(x-y) for field in ['error','position'] for x,y in zip(primary[field],ind[field]))<=1e-8,trace=max(traces.values())<=1e-9,resources=time.monotonic()-tick<=300 and full.memory.memory_bytes()<=2*1024**3)
                artifact=out/f'rational-a{amp:g}-n2-h{h:g}.npz';full.save(artifact,a,s);record=dict(amplitude=amp,n=2,h=h,primary=primary,independent=ind,traces=traces,rcond=s['rcond'],checks={k:bool(v) for k,v in checks.items()},seconds=time.monotonic()-tick,artifact=artifact.relative_to(ROOT).as_posix(),artifact_sha256=full.digest(artifact));records.append(record);print(record,flush=True);del a,s;gc.collect();assert all(checks.values())
            except Exception as exc:stops.append(dict(amplitude=amp,h=h,reason=repr(exc)));break
        if stops:break
    accuracy=[]
    if len(records)==4 and not stops:
        lookup={(r['amplitude'],r['h']):r for r in records}
        for amp in [0.,.2]:
            for field in ['error','position']:
                for target in range(3):
                    coarse=[lookup[amp,h]['independent'][field][target] for h in [.01,.0001]];fine=[next(r for r in base['records'] if r['mode']=='rational' and r['amplitude']==amp and r['h']==h)['independent'][field][target] for h in [.01,.0001]]
                    floor=max(coarse+fine)<1e-12;convergence=all(y<x for x,y in zip(coarse,fine));ratio=coarse[1]/max(coarse[0],1e-300)
                    accuracy.append(dict(amplitude=amp,field=field,target=target,n2=coarse,n4=fine,numerical_floor=floor,spatial_decrease=convergence,thin_ratio_n2=ratio,passed=(convergence or floor) and ratio<=2))
    paths=[parent,Path(__file__),ROOT/'rd/structure/main-shell/paired-geometry/records/paired-geometry-ancf-spatial-protocol.md',Path(full.__file__)]
    result=dict(passed=len(records)==4 and not stops and all(r['passed'] for r in accuracy) and time.monotonic()-start<=300,records=records,accuracy=accuracy,stops=stops,seconds=time.monotonic()-start,peak_bytes=full.memory.memory_bytes(),main_shell_adopted=False,sources={p.relative_to(ROOT).as_posix():full.digest(p) for p in paths});dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print('PASSED',result['passed'],flush=True)

if __name__=='__main__':main()
