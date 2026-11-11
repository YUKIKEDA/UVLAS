"""Frozen necessary improvement gates and retained v1/v2 field comparison."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json
import numpy as np
import moment_enriched_thin_v2 as t
ROOT=rd_paths.ROOT;D=rd_paths.RECORDS
def read(name):return json.loads((D/name).read_text())
def key(r):return r['kappa'],r['n'],r['amplitude']
def main():
    dest=D/'moment-enriched-thin-assessment.json';assert not dest.exists()
    primary=read('moment-enriched-thin-global-v2-result.json');old=read('reference-work-thin-global-v4-result.json');audit=read('moment-enriched-thin-global-v2-audit.json')
    comparisons=[]
    for r in primary['records']:
        b=next(v for v in old['records'] if key(v)==key(r));q,e=r['primary']['error'];oq,oe=b['primary']['error']
        comparisons.append(dict(kappa=r['kappa'],n=r['n'],amplitude=r['amplitude'],quartic=q,exponential=e,previous_quartic=oq,previous_exponential=oe,
            required=r['kappa']==.4 and r['amplitude']==.2,passed=(q<=1e-8 and e<=oe) if r['kappa']==.4 and r['amplitude']==.2 else True))
    original=read('moment-enriched-thin-global-result.json');first=original['records'][0];nextfirst=primary['records'][0]
    a=np.load(ROOT/first['artifact']);b=np.load(ROOT/nextfirst['artifact']);grad=jump=0.
    for ids,av,bv in zip(a['cells'],a['values'],b['values']):
        op=t.operators(a['ref'][ids].reshape(-1));s=op['state'];ar=av.reshape(3,s['nscalar'],3,2);br=bv.reshape(3,s['nscalar'],3,2)
        for k in range(3):
            ga=np.einsum('pjn,nak->pjak',s['Gwref'],ar[k]);gb=np.einsum('pjn,nak->pjak',s['Gwref'],br[k]);grad=max(grad,float(np.linalg.norm(ga-gb)/max(1.,np.linalg.norm(gb))))
            for face in s['faces']:
                ja=np.einsum('pn,nak->pak',face['jump'],ar[k]);jb=np.einsum('pn,nak->pak',face['jump'],br[k]);jump=max(jump,float(np.linalg.norm(ja-jb)/max(1.,np.linalg.norm(jb))))
    checks=dict(local=read('moment-enriched-thin-local-result.json')['all_checks_pass'],director=read('moment-enriched-director-image-result.json')['all_checks_pass'],global_numerics=primary['all_analysis_checks'],tensor_audit=audit['all_checks_pass'],
        comparisons=len(comparisons)==5 and all(x['passed'] for x in comparisons),basis_version_agreement=max(grad,jump)<=1e-8)
    ps=[Path(__file__),D/'moment-enriched-thin-global-v2-result.json',D/'reference-work-thin-global-v4-result.json',D/'moment-enriched-thin-global-v2-audit.json',D/'moment-enriched-thin-global-result.json',D/'moment-enriched-thin-local-result.json',D/'moment-enriched-director-image-result.json']
    dest.write_text(json.dumps(dict(comparisons=comparisons,retained_first_case_basis_comparison=dict(gradient=grad,stabilization=jump),checks=checks,all_necessary_checks_pass=all(checks.values()),finite_h_adoption=False,main_shell_adopted=False,geometry_conditional=True,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in ps}),indent=2),encoding='utf-8')
    print('assessment',checks,comparisons,'basis',grad,jump,flush=True)
if __name__=='__main__':main()
