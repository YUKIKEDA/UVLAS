from pathlib import Path
import sys,numpy as np
from diagnose_high_frequency_symmetric import unique_field
from run_high_frequency_final_bounded import bounded_schur
from uvlm_planar_decimal_residual import field as decimal_field,complex_dot
def physical_residual(g,tail,full,selected,eps,WW,rr):
    P=g['P'][selected];N=g['normal'][selected];known=complex_dot(WW,full[-g['n']:]);u=unique_field(g,tail,full,P,N,eps);R=u+known-rr;errors=np.linalg.norm(R,axis=0)/np.linalg.norm(rr,axis=0);print('UNCONDENSED_UNIQUE_FLOAT',errors.tolist(),flush=True)
    if max(errors)>=1e-11:
        u=decimal_field(g,tail,full,P,N,eps);R=u+known-rr;errors=np.linalg.norm(R,axis=0)/np.linalg.norm(rr,axis=0);print('UNCONDENSED_UNIQUE_DECIMAL70',errors.tolist(),flush=True)
    return errors
p=Path(__file__).with_name('run_high_frequency_final_symmetric.py');source=p.read_text(encoding='utf-8').replace("build/analysis/high-frequency-final'","build/analysis/high-frequency-final-qualified'").replace('assert max(sample)<1e-11;transfer=',"ring_sample=sample.copy();sample=physical_residual(g,tail,full,selected,eps,WW,rr);assert max(sample)<1e-11;transfer=").replace('uncondensed_sample_residuals=sample.tolist(),','uncondensed_sample_residuals=sample.tolist(),ring_column_diagnostic_residuals=ring_sample.tolist(),');namespace={'__name__':'qualified_experiment','__file__':str(p),'physical_residual':physical_residual};exec(compile(source,str(p),'exec'),namespace);namespace['schur_solve']=bounded_schur
if __name__=='__main__':namespace['main']()
