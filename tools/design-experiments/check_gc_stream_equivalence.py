"""Compare preserved non-streamed and streamed states in physical fields."""
from pathlib import Path
import hashlib,json
import numpy as np
import geometry_compatible_ancf_stream as f
ROOT=Path(__file__).resolve().parents[2]

def main():
    p=ROOT/'docs/verification/geometry-compatible-ancf-candidate-result.json';q=ROOT/'docs/verification/geometry-compatible-ancf-candidate-stream-result.json'
    old=json.loads(p.read_text());new=json.loads(q.read_text());select=lambda r:r['mode']=='standard' and r['amplitude']==.2 and r['n']==4 and r['h']==.01
    r=next(r for r in old['records'] if select(r));s=next(r for r in new['records'] if select(r));ap=ROOT/r['artifact'];bp=ROOT/s['artifact']
    with np.load(ap) as d:oa={k:d[k].copy() for k in d.files}
    with np.load(bp) as d:na={k:d[k].copy() for k in d.files}
    geometry_equal=all(np.array_equal(oa[k],na[k]) for k in oa if k not in ['coef','beta'])
    geo=f.c.prepare_geometry(.4,4,.2,13);geo.update(refs=list(oa['refs']),coeff=oa['geo_coeff'],V=oa['geo_V'],s=oa['geo_s'])
    for j,e in enumerate(geo['enrich']):
        for key in ['coeff','V','s']:e[key]=oa[f'enrich{j}_{key}']
    dc=na['coef']-oa['coef'];db=na['beta']-oa['beta'];ns=len(dc)//3
    values=[];strains=[];normv=[];norme=[]
    points=np.array([[x,y,z] for x in np.linspace(-.97,.97,7) for y in np.linspace(-.93,.93,7) for z in [-1.,-.41,0.,.37,1.]])
    for cell,ref in enumerate(geo['refs']):
        N,G,B,det,X=f.c.kinematics(4,cell,ref,.01,points,geo);E,U=f.bending(ref,points,.01,.4,4)
        values.append(((N@dc.reshape(ns,9)).reshape(len(points),3,3)+U@db).ravel());strains.append((B@dc+E@db).ravel())
        normv.append(((N@oa['coef'].reshape(ns,9)).reshape(len(points),3,3)+U@oa['beta']).ravel());norme.append((B@oa['coef']+E@oa['beta']).ravel())
    position=float(np.linalg.norm(np.concatenate(values))/np.linalg.norm(np.concatenate(normv)));strain=float(np.linalg.norm(np.concatenate(strains))/np.linalg.norm(np.concatenate(norme)))
    errors=max(abs(x-y) for field in ['error','position'] for x,y in zip(r['independent'][field],s['independent'][field]))
    paths=[p,q,ap,bp,Path(__file__)]
    result=dict(passed=geometry_equal and max(position,strain,errors)<=1e-8,geometry_equal=geometry_equal,relative_position_difference=position,relative_strain_difference=strain,error_difference=errors,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    dest=ROOT/'docs/verification/geometry-compatible-ancf-stream-equivalence.json';assert not dest.exists();dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(result)

if __name__=='__main__':main()
