import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import json,hashlib,numpy as np
ROOT=rd_paths.ROOT;OUT=ROOT/'build/analysis/unsteady-2-3'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    artifacts={};records=[];limits={};unc={};finite={}
    def get(label):
        p=OUT/(label+'.json');artifacts[p.relative_to(ROOT).as_posix()]=sha(p);snap=p.with_suffix('.npz');artifacts[snap.relative_to(ROOT).as_posix()]=sha(snap);return np.array([complex(*v) for v in json.loads(p.read_text())['transfer']])
    def extrap(a,b,c):
        r=abs(b-a)/abs(c-b);p=np.log2(r);direction=np.real((b-a)*np.conj(c-b))/(abs(b-a)*abs(c-b));assert np.all((p>.7)&(p<1.3)) and min(direction)>.99;return c+(c-b)/(r-1),p,r
    for side,old in [('candidate','candidate'),('reference','ptera-frequency')]:
        c8=get(f'{side}-continuum-m8-n96-q12');c16=get(f'{old}-m16-n96');c32=get(f'{side}-extension-m32-n96');c64=get(f'{side}-extension-m64-n96');n48=get(f'{old}-m16-n48');n192=get(f'{side}-continuum-m16-n192-q12');n384=get(f'{side}-extension-m16-n384');f=get(f'{side}-extension-m64-n192');x32=get(f'{side}-extension-m32-n192');oldfine=f.copy();quad=np.zeros(2)
        if side=='reference':f=get('reference-continuum-m64-n192-q16');quad=abs(f-oldfine)
        em0,pm0,_=extrap(c8,c16,c32);em1,pm1,rm=extrap(c16,c32,c64);en0,pn0,_=extrap(n48,c16,n192);en1,pn1,rn=extrap(c16,n192,n384);_,pmfine,rmfine=extrap(n192,x32,f);assert max(abs(pm0-pm1))<.25 and max(abs(pn0-pn1))<.25
        cm=(f-x32)/(rmfine-1);cn=(n384-n192)+(n384-n192)/(rn-1);limit=f+cm+cn;cross=abs((f-x32)-(c64-c32));estimate=(1.5*(abs(em1-em0)+abs(en1-en0))+cross+2*quad)/abs(limit);estimate+=3e-6+3e-6
        limits[side]=limit;unc[side]=estimate;finite[side]=f;records.append(dict(case='two_direction_continuum_estimate',side=side,chord_orders_coarse=pm0.tolist(),chord_orders_fine=pm1.tolist(),chord_order_at_n192=pmfine.tolist(),span_orders_coarse=pn0.tolist(),span_orders_fine=pn1.tolist(),limit=[[float(z.real),float(z.imag)] for z in limit],estimated_relative_uncertainty=estimate.tolist(),raw_finite_to_estimated_limit=abs((f-limit)/limit).tolist(),fine_quadrature_change=(quad/abs(limit)).tolist(),reference_space_target_met=bool(max(estimate)<.0025),rigorous_error_bound=False))
    difference=abs((limits['candidate']-limits['reference'])/limits['reference']);raw=abs((finite['candidate']-limits['reference'])/limits['reference']);total=raw+unc['reference'];records.append(dict(case='independent_continuum_and_raw_finite_comparison',estimated_limit_difference=difference.tolist(),candidate_raw_to_reference_limit=raw.tolist(),with_reference_uncertainty=total.tolist(),one_percent_with_reference_uncertainty=bool(max(total)<.01)))
    sources={p.relative_to(ROOT).as_posix():sha(p) for p in [Path(__file__),ROOT/'rd/aero/unsteady-reference/experiments/run_frequency_space_extension.py',ROOT/'rd/aero/unsteady-reference/experiments/run_frequency_continuum_checks.py',ROOT/'rd/aero/unsteady-reference/records/unsteady-space-extension-protocol.md',ROOT/'rd/aero/unsteady-reference/records/unsteady-continuum-estimate-protocol.md',ROOT/'rd/aero/unsteady-reference/records/unsteady-initial-campaign-result.json']};reference_qualified=bool(max(unc['reference'])<.0025);r=dict(passed=reference_qualified,scope='k=0.2 infinitesimal heave FW-V01; continuum reference estimate',reference_continuum_qualified=reference_qualified,raw_candidate_original_space_allocation_met=False,uvlm_2_3_complete=False,records=records,sources=sources,artifacts=artifacts);dest=ROOT/'rd/aero/unsteady-reference/records/unsteady-continuum-estimate-result.json';assert not dest.exists();dest.write_text(json.dumps(r,indent=2),encoding='utf-8');print(json.dumps({k:v for k,v in r.items() if k not in ['sources','artifacts']},indent=2))
if __name__=='__main__':main()
