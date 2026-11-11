"""Audit original AVL residuals, segment forces and instrumentation effects."""
from pathlib import Path
import hashlib
import json
import math
import re
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
FOLDER=ROOT/'docs/verification'
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def fields(path):
    values={}
    for line in path.read_text().splitlines():
        if '|' not in line:continue
        left,right=line.split('|',1)
        try:values[right.strip()]=[float(v) for v in left.split()]
        except ValueError:pass
    return values

def coefficients(v):
    return dict(force=[v["CXtot, Cltot, Cl'tot"][0],v['CYtot, Cmtot'][0],v["CZtot, Cntot, Cn'tot"][0]],
        moment=[v["CXtot, Cltot, Cl'tot"][1],v['CYtot, Cmtot'][1],v["CZtot, Cntot, Cn'tot"][1]],
        lift=v['CLtot'][0],drag=v['CDtot'][0])

def main():
    parent_path=FOLDER/'avl-stable-reference-runs.json'
    parent=json.loads(parent_path.read_text());assert parent['completed']
    records=[];by_key={}
    for run in parent['runs']:
        folder=ROOT/run['directory']
        total=fields(folder/'total.txt');coefs=coefficients(total)
        log=(folder/'console.txt').read_text()
        flags={name:re.findall(re.escape(label)+r'\s*:\s*([TF])',log)[-1] for name,label in
            [('trailing','T rail.leg forces'),('viscous','V iscous forces'),('body','B ody forces')]}
        item=dict(variant=run['variant'],case=run['case'],m=run['m'],n=run['n'],alpha=run['alpha'],
            coefficients=coefs,options=flags,
            measurement_recompute_unchanged=all(digest(folder/f'before-{name}.txt')==digest(folder/f'{name}.txt') for name in ['total','strips','elements']))
        if run['variant']!='baseline':
            rows=np.loadtxt(folder/'uv-rows.txt',skiprows=2)
            forces=np.loadtxt(folder/'uv-forces.txt',skiprows=2)
            matrix_data=np.loadtxt(folder/'uv-matrix.txt',skiprows=1)
            n=run['m']*run['n'];assert rows.shape==(n,21)
            assert matrix_data.shape==(n*n,3) and forces.shape==(3*n,9)
            assert np.array_equal(rows[:,0],np.arange(1,n+1))
            assert np.array_equal(matrix_data[:,:2],np.array([(i,j) for i in range(1,n+1) for j in range(1,n+1)]))
            matrix=matrix_data[:,2].reshape(n,n);gamma=rows[:,13]
            known=rows[:,14:17];normal=rows[:,4:7]
            rhs=-np.einsum('ij,ij->i',normal,known)
            residual=np.array([math.fsum(float(matrix[i,j])*float(gamma[j]) for j in range(n))-rhs[i] for i in range(n)])
            native=rows[:,20]
            scale=max(np.linalg.norm(known,axis=1).max(),1e-300)
            f=np.array([math.fsum(forces[:,6+j]) for j in range(3)])
            origin=np.array(total['Xref, Yref, Zref'])
            torques=np.cross(forces[:,3:6]-origin,forces[:,6:9])
            torque=np.array([math.fsum(torques[:,j]) for j in range(3)])
            s,c,b=total['Sref, Cref, Bref']
            expected_m=torque/(s*np.array([b,c,b]));expected_f=f/s
            tuples={(int(row[0]),int(row[2])) for row in forces}
            item.update(dict(residual_relative=float(np.max(abs(residual))/scale),
                native_residual_relative=float(np.max(abs(native))/scale),
                residual_difference=float(np.max(abs(residual-native))),
                force_difference=float(np.max(abs(expected_f-coefs['force']))),
                moment_difference=float(np.max(abs(expected_m-coefs['moment']))),
                force_records=len(forces),force_ids_complete=tuples=={(i,k) for i in range(1,n+1) for k in range(3)},
                all_finite=bool(all(np.all(np.isfinite(x)) for x in [rows,forces,matrix_data])),
                gamma_over_U_l2=float(np.linalg.norm(gamma)),coefficient_force_from_records=expected_f.tolist(),
                coefficient_moment_from_records=expected_m.tolist()))
        records.append(item)
        by_key[(run['variant'],run['case'],run['m'],run['n'],run['alpha'])]=(run,item)
    unchanged=[]
    for case in ['FW-V01','FW-V02']:
        for m,n in [(4,12),(8,24)]:
            for alpha in [-2,-1,0,1,2]:
                raw,_=by_key[('baseline',case,m,n,alpha)]
                changed,_=by_key[('measured',case,m,n,alpha)]
                same=all(digest(ROOT/raw['directory']/f'{name}.txt')==digest(ROOT/changed['directory']/f'{name}.txt') for name in ['total','strips','elements'])
                unchanged.append(dict(case=case,m=m,n=n,alpha=alpha,outputs_equal=same))
    official=json.loads((FOLDER/'avl-reference-execution-corrected-result.json').read_text())
    official_cases=[]
    for case,label in [('FW-V01','planar'),('FW-V02','dihedral')]:
        old=next(r for r in official['runs'] if r['label']==label)
        old_values=coefficients(fields(ROOT/old['directory']/'total.txt'))
        _,current=by_key[('baseline',case,4,12,1)]
        diff=max(abs(x-y) for key in ['force','moment'] for x,y in zip(old_values[key],current['coefficients'][key]))
        official_cases.append(dict(case=case,maximum_absolute_difference=diff,eligible=diff<=1e-7))
    sensitivity=[]
    for case in ['FW-V01','FW-V02']:
        for m,n in [(4,12),(8,24)]:
            for delta in [1,2]:
                slopes={}
                for variant in ['measured','core0']:
                    plus=by_key[(variant,case,m,n,delta)][1]['coefficients']
                    minus=by_key[(variant,case,m,n,-delta)][1]['coefficients']
                    denom=2*np.deg2rad(delta)
                    slopes[variant]=dict(lift=(plus['lift']-minus['lift'])/denom,
                        moment_y=(plus['moment'][1]-minus['moment'][1])/denom)
                lift=abs(slopes['measured']['lift']-slopes['core0']['lift'])/abs(slopes['measured']['lift'])
                moment=abs(slopes['measured']['moment_y']-slopes['core0']['moment_y'])/max(abs(slopes['measured']['moment_y']),.1)
                sensitivity.append(dict(case=case,m=m,n=n,delta_deg=delta,slopes=slopes,lift_difference=lift,moment_difference=moment,eligible=max(lift,moment)<=.0005))
    measured=[r for r in records if r['variant']!='baseline']
    checks=dict(all_runs_complete=len(records)==60,
        options=all(r['options']==dict(trailing='T',viscous='F',body='F') for r in records),
        measurement_recompute=all(r['measurement_recompute_unchanged'] for r in records),
        original_vs_measured=all(r['outputs_equal'] for r in unchanged),
        original_vs_official=all(r['eligible'] for r in official_cases),
        residual=all(max(r['residual_relative'],r['native_residual_relative'])<=1e-11 for r in measured),
        residual_agreement=all(r['residual_difference']<=1e-12 for r in measured),
        force_sum=all(r['force_difference']<=1e-10 for r in measured),
        moment_sum=all(r['moment_difference']<=1e-10 for r in measured),
        records_complete=all(r['force_ids_complete'] and r['all_finite'] for r in measured),
        core_sensitivity=all(r['eligible'] for r in sensitivity),
        runtime_budget=parent['seconds']<=300,
        memory_budget=all(r['peak_resident_bytes']<2*1024**3 for r in parent['runs']),
        output_budget=sum(v['bytes'] for r in parent['runs'] for v in r['files'].values())<=256*1024**2)
    result=dict(protocol='docs/verification/avl-source-measurement-experiment.md',source_sha256=digest(Path(__file__)),
        parent_result_sha256=digest(parent_path),official_result_sha256=digest(FOLDER/'avl-reference-execution-corrected-result.json'),
        cases=records,output_comparisons=unchanged,official_comparisons=official_cases,core_sensitivity=sensitivity,
        checks=checks,all_criteria_pass=all(checks.values()),
        limitations='reference instrumentation only; no fixed-strip conservation, full grid convergence or UVLAS load comparison')
    (FOLDER/'avl-stable-reference-result.json').write_text(json.dumps(result,indent=2,allow_nan=False,default=lambda v:v.item() if isinstance(v,np.generic) else str(v))+'\n')
    print(json.dumps(dict(checks=checks,max_residual=max(r['residual_relative'] for r in measured),
        max_force_difference=max(r['force_difference'] for r in measured),
        max_moment_difference=max(r['moment_difference'] for r in measured),
        max_core_sensitivity=max(max(r['lift_difference'],r['moment_difference']) for r in sensitivity))))
    if not all(checks.values()):raise SystemExit(1)

if __name__=='__main__':main()
