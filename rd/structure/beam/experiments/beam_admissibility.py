"""Conservative rigid-section tube and explicit unsupported-warping guards."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
import beam_se3_v2 as beam
import beam_section_input as inputs
def prepare(request):
    section=inputs.prepare(request);bounds=np.asarray(request.get('section_envelope'),float)
    if bounds.shape!=(2,) or not np.isfinite(bounds).all() or np.any(bounds<=0):raise ValueError('missing_or_invalid_section_envelope')
    section['section_envelope']=bounds.copy();section['section_envelope'].flags.writeable=False
    return section
def regularity(r,R,L,section):
    beam.validate(r,R,L);d=beam.float_array(beam.relative(r[0],R[0],r[1],R[1]))/L
    y,z=section['section_envelope'];lower=float(d[0]-abs(d[4])*z-abs(d[5])*y)
    if lower<=0:raise ValueError('section_tube_regularity_uncertified')
    return dict(jacobian_lower=lower,sufficient_not_necessary=True,formula='Gamma_x+kappa_y*z-kappa_z*y')
