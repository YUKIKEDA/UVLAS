"""Equivalent block equilibration of the 174 basic plus 3 completion equations."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
from scipy.linalg import qr,svd,solve_triangular
import moment_enriched_thin_v2 as t
def nullspace(A,independent=False):
    singular=svd(A,compute_uv=False,check_finite=False);raw=t.rank_info(singular)
    if len(set(raw))==1:
        Z,info=t.nullspace(A,independent=independent);info.update(raw_ranks=raw,block_equilibrated=False);return Z,info
    if A.shape[0]!=177:raise ValueError('MMS1_unexpected_ambiguous_constraint_size:'+str(A.shape))
    _,_,pivot=qr(A,mode='economic',pivoting=True,check_finite=False);strong=pivot[:174];free=pivot[174:];G=A[:,strong];sv=svd(G,compute_uv=False,check_finite=False);sr=t.rank_info(sv)
    if sr!=[174]*3:raise ValueError('MMS1_strong_block_rank:'+str(sr))
    Q,R=qr(G,mode='full',check_finite=False);top=Q[:,:174].T@A[:,free];weak=Q[:,174:].T@A[:,free];scale=np.linalg.norm(weak)
    if scale==0:raise ValueError('MMS1_zero_completion_block')
    Zfree,wi=t.nullspace(weak/scale,independent=independent);Z=np.zeros((A.shape[1],Zfree.shape[1]));Z[free]=Zfree;Z[strong]=-solve_triangular(R[:174],top@Zfree,check_finite=False);Z=qr(Z,mode='economic',check_finite=False)[0]
    residual=float(np.linalg.norm(A@Z)/max(1.,np.linalg.norm(A)));weak_res=float(np.linalg.norm(weak@Z[free])/scale)
    if residual>1e-9 or weak_res>1e-9:raise ValueError('MMS1_equivalent_block_residual')
    rank=174+wi['ranks'][0]
    return Z,dict(raw_ranks=raw,ranks=[rank]*3,columns=A.shape[1],nullity=Z.shape[1],residual=residual,block_equilibrated=True,strong_ranks=sr,strong_rcond=float(sv[-1]/sv[0]),completion=wi,completion_scale=float(scale),completion_residual=weak_res)
