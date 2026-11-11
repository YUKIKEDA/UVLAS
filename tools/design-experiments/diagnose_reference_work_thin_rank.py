"""Reproduce the rank sensitivity that motivated symbolic compatibility."""
from pathlib import Path
import hashlib,json
import numpy as np
from scipy.linalg import svd
import reference_work_thin_limit as t
ROOT=Path(__file__).resolve().parents[2]
def main():
    p=ROOT/'docs/verification/reference-work-thin-rank-diagnostic.json';assert not p.exists();g=t.m.mesh(.4,2,.2);q=g['ref'][g['cells'][0]].reshape(-1);op=t.operators(q);C=t.constraints(op);scale=1/np.linalg.norm(C,axis=0);R=t.f.qrtools.rfactor(np.asfortranarray(C*scale));s=svd(R,compute_uv=False)
    p.write_text(json.dumps(dict(monolithic_cutoff_ranks=t.rank_info(s),relative_cutoffs=[1e-9,1e-10,1e-11],singular_values_900_to_908=s[900:909].tolist(),largest_singular_value=float(s[0]),
        chosen_adoption_rank=None,exact_rank_of_rounded_geometry_certified=False,analytic_geometry_certificate_is_a_different_statement=True,
        sources={Path(__file__).relative_to(ROOT).as_posix():hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}),indent=2),encoding='utf-8')
if __name__=='__main__':main()
