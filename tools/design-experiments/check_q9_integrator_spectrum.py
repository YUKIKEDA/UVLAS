"""Supplementary finite-frequency analysis, protocol section 6; no nonlinear rerun."""
import hashlib
import json
from pathlib import Path
import numpy as np
import check_q9_integrator_candidates as candidates

ROOT=Path(__file__).resolve().parents[2]
source=Path(__file__)
rows=[]
for rho in [1.,.8,.5,0.]:
    for eta in [.01,.1,1.,10.]:
        values=np.linalg.eigvals(candidates.amplification(rho,eta))
        oscillatory=values[values.imag>1e-12]
        row=dict(rho=rho,h_omega=eta,spectral_radius=float(np.abs(values).max()),
            eigenvalues=[[float(v.real),float(v.imag)] for v in values])
        if len(oscillatory)==1:
            lam=oscillatory[0]
            row.update(oscillatory_amplitude_per_step=float(abs(lam)),frequency_ratio=float(np.angle(lam)/eta))
        else:row['oscillatory_branch']='absent_or_ambiguous'
        rows.append(row)
input_path=ROOT/'docs/verification/q9-integrator-candidate-result.json'
previous=json.loads(input_path.read_text(encoding='utf-8'))
orders=[]
for rho in [.8,.5,0.]:
    cases={c['n']:c for c in previous['nonlinear_cases'] if c['rho']==rho}
    orders.append(dict(rho=rho,observed_orders_640_1280={k:float(np.log2(cases[640]['errors'][k]/cases[1280]['errors'][k]))
        for k in ['displacement_mass_max','velocity_scaled_mass_max','physical_acceleration_relative','pk2_increment_relative','low4_phase_space_max']}))
result=dict(protocol='docs/verification/q9-integrator-candidate-comparison.md#6',source_sha256=candidates.sha(source),
    dependency_sha256=candidates.sha(source.with_name('check_q9_integrator_candidates.py')),
    nonlinear_result_sha256=candidates.sha(input_path),finite_frequency_cases=rows,nonlinear_observed_orders=orders)
(ROOT/'docs/verification/q9-integrator-spectrum-result.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps(dict(cases=[r for r in rows if r['h_omega'] in [.1,1.]],orders=orders)))
