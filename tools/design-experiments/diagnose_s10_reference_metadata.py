from pathlib import Path
import json,hashlib,zipfile,io,re
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'build/analysis/s10-reassessment'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def real(x):
    x=x.strip().replace('D','E').replace('d','e')
    if not x:return 0.
    if 'e' not in x.lower():x=re.sub(r'(?<=[0-9.])([+-][0-9]+)$',r'E\1',x)
    return float(x)
def mass_cards(path):
    cards=[];current=None
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith('$'):continue
        name=line[:8].strip()
        if name.startswith('+') or name=='*':
            if current is not None:current[1].extend([line[i:i+8] for i in range(8,72,8)])
        else:
            current=[name,[line[i:i+8] for i in range(8,72,8)]];cards.append(current)
    values=[]
    for name,f in cards:
        if name=='CONM2':values.append(dict(eid=int(f[0]),node=int(f[1]),cid=int(f[2] or 0),mass=real(f[3]),offset=[real(s) for s in f[4:7]]))
    return dict(count=len(values),sum_mass=sum(x['mass'] for x in values),largest=sorted(values,key=lambda x:x['mass'],reverse=True)[:8],scope='CONM2 only; not total FE model mass; read-only metadata diagnostic, not an importer')
def main():
    dest=ROOT/'docs/verification/s10-reference-metadata-result.json';assert not dest.exists();sources={};artifacts={};record={}
    with zipfile.ZipFile(ROOT/'build/reference-audit/s10-analyzed.zip') as z:
        record['published_modes']={}
        for config in ['LE','TE']:
            d=loadmat(io.BytesIO(z.read(f'Analyzed_data/modes_{config}.mat')),simplify_cells=True)[f'modes_{config}'];freq=np.asarray(d['natural_freq'],float);record['published_modes'][config]=dict(frequencies=freq.tolist(),marker_coordinates_shape=list(d['marker_position_in_undeformed_wing_xyz'].shape),mode_count=len(d['natural_modes_Z']))
            assert len(freq)==3 and np.isfinite(freq).all()
    static=next(OUT.glob('*4012*mat'));d=loadmat(static,simplify_cells=True)['raw_data_organized'];q=np.asarray(d['DynP_1000hz']);v=np.asarray(d['vel_1000hz']);density=2*q*98100/v**2;column_error=float(np.max(np.abs(q-d['tunnel_data1000hz'][:,10])));assert column_error==0.
    record['dynamic_pressure']=dict(run=static.name,column_11_exact_match=True,factor_to_Pa=98100.,implied_density_quantiles=np.quantile(density,[0,.5,1]).tolist(),meaning='internal unit consistency only; not calibration uncertainty or independently measured density')
    record['model_conm2']={config:mass_cards(ROOT/f'build/reference-audit/S10_{config}.nas') for config in ['LE','TE']}
    missing={}
    for name in ['strain_sensor_x_location_LE','strain_sensor_x_location_TE']:missing[name]=name not in d
    inventory=json.loads((OUT/'public-inspection.json').read_text());record['sensor_location_availability']=dict(missing_in_inspected_static_fields=missing,sensor_loc_m_in_public_raw_zip=any(Path(x['name']).name.lower()=='sensor_loc.m' for x in inventory['raw_inventory']),meaning='availability audit; do not invent x from y')
    for name in ['Run_Numbers.xlsx','public-inspection.json','case-input-audit.json']+[p.name for p in OUT.glob('*.mat')]:p=OUT/name;artifacts[p.relative_to(ROOT).as_posix()]=sha(p)
    for name in ['build/reference-audit/s10-analyzed.zip','build/reference-audit/s10-metadata.json','build/reference-audit/s10-readme.pdf','build/reference-audit/s10-readme.txt','build/reference-audit/S10_LE.nas','build/reference-audit/S10_TE.nas']:artifacts[name]=sha(ROOT/name)
    for name in ['inspect_s10_public_resources.py','reassess_s10_case_inputs.py',Path(__file__).name]:sources['tools/design-experiments/'+name]=sha(ROOT/'tools/design-experiments'/name)
    result=dict(passed=True,experimental_reference_qualified=False,product_acceptance_complete=False,scope='public metadata, units and availability; no experimental accuracy pass',records=record,remaining=['measurement/calibration/time/geometry uncertainty','camera marker and coordinate registration','configuration mass/fuse correspondence','strain sensor x and sensor axes','pre-flutter stationary small-amplitude intervals and attached-flow applicability'],sources=sources,artifacts=artifacts);dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(record,indent=2))
if __name__=='__main__':main()
