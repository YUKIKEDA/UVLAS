"""Explicit sectional constitutive/mass data; no restrained-warping inference."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
import beam_se3_v2 as beam
ORDER='gamma_x_gamma_y_gamma_z_kappa_x_kappa_y_kappa_z'
def prepare(request):
    if request.get('units')!='SI' or request.get('order')!=ORDER:raise ValueError('units_or_order')
    if request.get('warping')!='saint_venant_free' or request.get('restrained_warping',False):raise ValueError('unsupported_restrained_warping')
    if not request.get('source_id'):raise ValueError('missing_section_provenance')
    L=float(request['reference_scale']);K=np.asarray(request['K'],float);beam.constitutive(K,L)
    e0=np.asarray(request['initial_generalized_strain'],float)
    if e0.shape!=(6,) or not np.isfinite(e0).all():raise ValueError('invalid_initial_generalized_strain')
    mass=beam.spatial_mass(float(request['mass_per_length']),np.array(request['centroid'],float),np.array(request['centroid_inertia_per_length'],float))
    out=dict(K=K.copy(),mass=mass,e0=e0.copy(),source_id=request['source_id'],warping=request['warping'],local_stress_recovery=request.get('local_stress_recovery','unavailable'))
    for name in ['K','mass','e0']:out[name].flags.writeable=False
    return out
