"""Negative controls only: selectively remove r/d intercell continuity."""
from pathlib import Path
import hashlib,json,itertools,time
import numpy as np
import check_q4_continuity_v2 as checker
import tying_mid_projection_q4_v2 as m
ROOT=Path(__file__).resolve().parents[2]

def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def relax(original,broken):
    ref=[];xy=[];field_ids=[];free=[];connect=[];offset=0
    for f in range(3):
        old_ids=original['field_ids'][f]
        if f in broken:
            ids=np.concatenate([cell[m.SLICES[f]] for cell in original['cells']])
            local=np.array([[offset+25*c+i for i in range(25)] for c in range(4)])
        else:
            ids=old_ids
            lookup={int(node):offset+i for i,node in enumerate(ids)}
            local=np.array([[lookup[int(node)] for node in cell[m.SLICES[f]]] for cell in original['cells']])
        coords=original['xy'][f][ids-old_ids[0]]
        ref.append(original['ref'][ids]);xy.append(coords);field_ids.append(np.arange(offset,offset+len(ids)))
        free.append(field_ids[-1][np.abs(coords[:,0])>1e-12]);connect.append(local);offset+=len(ids)
    gauges=[]
    if 0 in broken:
        for ids in connect[0]:
            if not np.any(np.abs(xy[0][ids,0])<1e-12):
                gauges.append(ids);free[0]=free[0][free[0]!=ids[0]]
    return dict(ref=np.concatenate(ref),xy=xy,field_ids=field_ids,free=free,r_translation_gauges=gauges,
                cells=np.array([np.concatenate([connect[f][c] for f in range(3)]) for c in range(4)]))

def main():
    dest=ROOT/'docs/verification/q4-continuity-v2-result.json';assert not dest.exists()
    out=ROOT/'build/analysis'/('q4-continuity-v2-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir()
    start=time.monotonic();records=[];stop=None
    try:
        for broken in [[],[0],[1],[0,1]]:
            mesh=relax(m.mesh(.4,2,.2),broken)
            for method in ['raw','conservative']:
                target=out/('broken-'+''.join(str(i) for i in broken)+'-'+method);target.mkdir()
                r=checker.shared(mesh,.4,.2,method,target)
                r['broken_fields']=broken;r['negative_control']=bool(broken)
                r['checks']['resources']=time.monotonic()-start<900
                records.append(r)
                print('broken',broken,method,'null',r['nullity'],{k:v['error'] for k,v in r['solutions'].items()},r['checks'],flush=True)
                if not r['checks']['resources']:raise RuntimeError('resource_gate_failed')
    except Exception as exc:stop=repr(exc)
    sources=[Path(__file__),Path(m.__file__),Path(checker.__file__),ROOT/'docs/verification/q4-continuity-v2-protocol.md',ROOT/'docs/verification/q4-continuity-v2-preparation.json']
    dest.write_text(json.dumps(dict(records=records,stop_reason=stop,physical_adoption=False,
        all_numerical_checks=len(records)==8 and all(all(r['checks'].values()) for r in records) and stop is None,
        sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start),indent=2),encoding='utf-8')

if __name__=='__main__':main()
