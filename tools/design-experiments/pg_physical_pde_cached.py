"""Reference-only bounded caches; preserve physical work and preparation costs."""
from pathlib import Path
from collections import OrderedDict
import hashlib,json,time
import numpy as np
import pg_physical_pde_constrained as route
pde=route.pde;ROOT=Path(__file__).resolve().parents[2]
_components=route.components;_basis=pde.loads.global_basis;_prepare=route.prepare
_values=OrderedDict();_basis_cache=OrderedDict();_value_bytes=0;_basis_bytes=0

def digest(a):return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()

def basis(x,k,order=32):
    global _basis_bytes
    x=np.asarray(x);unique,back=np.unique(x,return_inverse=True);key=(float(k),order,digest(unique))
    if key not in _basis_cache:
        value=_basis(unique,k,order)
        for a in value:a.setflags(write=False)
        size=sum(a.nbytes for a in value)
        while _basis_cache and _basis_bytes+size>64*1024**2:
            _,old=_basis_cache.popitem(last=False);_basis_bytes-=sum(a.nbytes for a in old)
        _basis_cache[key]=value;_basis_bytes+=size
    _basis_cache.move_to_end(key);K,Kx,U=_basis_cache[key];return K[back],Kx[back],U[:,back]

def components(a,cell,q,p):
    global _value_bytes
    if '_component_generation' not in a:
        g=a['geo'];e=g['enrich'][0];a['_component_generation']=tuple(digest(v) for v in [g['coeff'],g['V'],g['s'],e['coeff'],e['V'],e['s']])
    key=(a['n'],a['ns'],a['_component_generation'],digest(q),digest(p))
    if key not in _values:
        value=_components(a,cell,q,p);value.setflags(write=False)
        while _values and _value_bytes+value.nbytes>512*1024**2:
            _,old=_values.popitem(last=False);_value_bytes-=old.nbytes
        _values[key]=value;_value_bytes+=value.nbytes
    _values.move_to_end(key);return _values[key]

def prepare(k,n,amp,h,order=13):
    paths=[Path(__file__),Path(route.__file__),Path(pde.__file__),ROOT/'tools/design-experiments/paired_geometry_ancf.py',ROOT/'tools/design-experiments/paired_curvature_profiles.py',ROOT/'docs/verification/paired-geometry-ancf-candidate-audit.json']
    sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths};identity=dict(k=k,n=n,amplitude=amp,h=h,order=order,sources=sources);key=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest();folder=ROOT/'build/analysis/pg-physical-preparations';folder.mkdir(exist_ok=True);path=folder/(key+'.npz');meta=folder/(key+'.json');tick=time.monotonic()
    if path.exists():
        record=json.loads(meta.read_text());assert hashlib.sha256(path.read_bytes()).hexdigest()==record['artifact_sha256']
        with np.load(path) as d:
            a={x:d[x].copy() for x in ['Zs','Z','scale','R']};nr,nd,nc,ns=map(int,d['sizes']);geo=dict(refs=list(d['refs']),coeff=d['geo_coeff'].copy(),V=d['geo_V'].copy(),s=d['geo_s'].copy(),rank=len(d['geo_s']),metric=d['metric'].tolist(),n=n,enrich=[dict(coeff=d['enrich_coeff'].copy(),V=d['enrich_V'].copy(),s=d['enrich_s'].copy())]);a.update(k=k,n=n,amp=amp,h=h,geo=geo,nr=nr,nd=nd,nc=nc,ns=ns,rcond=record['rcond'],support_ranks=record['support_ranks'],commutation_error=record['commutation_error'],load_value_error=record['load_value_error'],prepare_seconds=record['prepare_seconds']+time.monotonic()-tick)
    else:
        a=_prepare(k,n,amp,h,order);g=a['geo'];e=g['enrich'][0];np.savez_compressed(path,Zs=a['Zs'],Z=a['Z'],scale=a['scale'],R=a['R'],sizes=[a['nr'],a['nd'],a['nc'],a['ns']],refs=np.array(g['refs']),geo_coeff=g['coeff'],geo_V=g['V'],geo_s=g['s'],metric=g['metric'],enrich_coeff=e['coeff'],enrich_V=e['V'],enrich_s=e['s']);a['prepare_seconds']=time.monotonic()-tick;record=dict(prepare_seconds=a['prepare_seconds'],rcond=a['rcond'],support_ranks=a['support_ranks'],commutation_error=a['commutation_error'],load_value_error=a['load_value_error'],artifact=path.relative_to(ROOT).as_posix(),artifact_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),sources=sources);meta.write_text(json.dumps(record,indent=2),encoding='utf-8')
    x=np.array([0.,.17,.43,.79,1.,.17]);actual=basis(x,k);expected=_basis(x,k);a['basis_cache_error']=max(float(np.linalg.norm(v-w)/max(1.,np.linalg.norm(w))) for v,w in zip(actual,expected));assert a['basis_cache_error']<=1e-12;a['preparation_manifest']=meta.relative_to(ROOT).as_posix();return a

pde.prepare=prepare;pde.components=components;pde.loads.global_basis=basis
