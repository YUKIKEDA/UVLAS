"""Saved GSP polynomial zero/first maps and physically weighted shared space."""
from pathlib import Path
import hashlib,json,itertools
import numpy as np
from scipy.linalg import solve
import gsp_q3_v2 as candidate
import q3_thin_expansion as common
from mq3i1_shared_space import solve_semidefinite_ball
ROOT=Path(__file__).resolve().parents[2]
KEYS=['q','frame_scale','length','center','T','active','column_norms','coeff','denominator']
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def key(q):return hashlib.sha256(q.tobytes()).hexdigest()
def rel(a,b):return float(np.linalg.norm(a)/max(1.,np.linalg.norm(b)))

def load_cache():
    cache={}
    for label in ['gsp-q3-local-v2','gsp-q3-additional-local']:
        path=ROOT/f'docs/verification/{label}-result.json';data=json.loads(path.read_text())
        audit=json.loads((ROOT/f'docs/verification/{label}-audit.json').read_text())
        assert data['all_checks_pass'] and audit['all_reproduction_checks'] and audit['primary_complete']
        for source,sha in data['sources'].items():assert digest(ROOT/source)==sha
        for item in data['states']:
            assert item['passed'];p=ROOT/item['artifact'];assert digest(p)==item['artifact_sha256']
            saved=np.load(p);state={k:saved[k] for k in KEYS}
            cache[key(state['q'])]=dict(state=state,a0=saved['a0'],ad=saved['ad'],artifact=item['artifact'],
                artifact_sha256=item['artifact_sha256'],local_result=path.relative_to(ROOT).as_posix())
    return cache

def validate(item):
    state=item['state'];p,w=candidate.BASE.quadrature(9)
    _,E0=candidate.linear_raw(state['q'],0.,p);_,Ec=candidate.linear_raw(state['q'],1e-25j,p)
    P0=candidate.project(state,0.,E0);Pc=candidate.project(state,1e-25j,Ec)
    ad=solve(P0['H'],Pc['rhs'].imag/1e-25-Pc['H'].imag/1e-25@P0['a'],assume_a='sym',check_finite=False)
    metrics=dict(zero_coefficients=rel(item['a0']-P0['a'],P0['a']),first_coefficients=rel(item['ad']-ad,ad))
    plane=np.array(list(itertools.product(*[np.polynomial.legendre.leggauss(11)[0]]*2)))
    fields=[]
    for z in [-1.,0.,1.]:
        pt=np.column_stack((plane,np.full(len(plane),z)));S=candidate.stress_field(state,pt)
        fields.append((np.einsum('ij,pjn->pin',candidate.CI,S@item['a0']),np.einsum('ij,pjn->pin',candidate.CI,S@item['ad'])))
    metrics['zero_thickness_dependence']=max(rel(fields[i][0]-fields[1][0],fields[1][0]) for i in [0,2])
    metrics['first_odd']=rel(fields[0][1]+fields[2][1],fields[2][1])
    metrics['first_zero']=rel(fields[1][1],fields[2][1])
    pp,ww=candidate.BASE.quadrature(11);S=candidate.stress_field(state,pp)
    R0=np.einsum('ij,pjn->pin',candidate.CI,S@item['a0']);Rd=np.einsum('ij,pjn->pin',candidate.CI,S@item['ad'])
    _,raw0=candidate.linear_raw(state['q'],0.,pp);_,rawc=candidate.linear_raw(state['q'],1e-25j,pp)
    mu,_=candidate.measure(state,0.,pp,ww);mc,_=candidate.measure(state,1e-25j,pp,ww)
    metrics['work_zero']=rel(np.einsum('p,pkn->kn',mu,R0-raw0),np.einsum('p,pkn->kn',mu,raw0))
    work1=np.einsum('p,pkn->kn',mu,Rd-rawc.imag/1e-25)+np.einsum('p,pkn->kn',mc.imag/1e-25,R0-raw0)
    metrics['work_first']=rel(work1,np.einsum('p,pkn->kn',mu,rawc.imag/1e-25))
    return metrics,{k:bool(v<=1e-9) for k,v in metrics.items()}

def extract(item,plane):
    state=item['state'];p0=np.column_stack((plane,np.zeros(len(plane))));p1=p0.copy();p1[:,2]=1.
    S0=candidate.stress_field(state,p0);S1=candidate.stress_field(state,p1)
    zero=np.einsum('ij,pjn->pin',candidate.CI,S0@item['a0'])
    first=2*np.einsum('ij,pjn->pin',candidate.CI,S1@item['ad'])
    _,J,_=candidate.kinematics(state['q'],0.,p0)
    return J,zero[:,:,:144],first[:,:,:144]

def assemble(k,n,amp,order,cache,maps):
    xy,cells,coords,_=common.element.mesh(k,n,amp);free=np.flatnonzero(np.arange(len(xy))%(3*n+1));lookup={int(v):i for i,v in enumerate(free)}
    x,w=np.polynomial.legendre.leggauss(order);plane=np.array(list(itertools.product(x,x)));weights=np.array([a*b for a,b in itertools.product(w,w)])
    rows=6*len(plane)*len(cells);A0=np.zeros((rows,6*len(free)));A1=np.zeros_like(A0);D1=np.zeros((rows,3*len(free)))
    targets={kind:np.zeros(rows) for kind in ['quadratic','quartic']};keys=[];zero_c=0.
    for cell,(_,_,ids) in enumerate(cells):
        q=coords[ids].reshape(144);name=key(q)
        if name not in cache:raise ValueError('unvalidated_reference_geometry:'+name)
        item=cache[name]
        if name not in maps:
            metrics,checks=validate(item)
            maps[name]={k:item[k] for k in ['artifact','artifact_sha256','local_result']}
            maps[name].update(metrics=metrics,checks=checks)
        if not all(maps[name]['checks'].values()):raise ValueError('zero_first_validation_failed:'+name)
        keys.append(name);J,B0,B1=extract(item,plane);weight=np.sqrt(weights*np.linalg.det(J))
        W0=np.einsum('ij,pjn->pin',candidate.CH,B0);W1=np.einsum('ij,pjn->pin',candidate.CH,B1)
        sl=slice(6*len(plane)*cell,6*len(plane)*(cell+1));zero_c=max(zero_c,float(np.linalg.norm(B0[:,:,common.CC])))
        for local,node in enumerate(ids):
            if int(node) not in lookup:continue
            index=lookup[int(node)]
            A0[sl,6*index:6*index+6]=(W0[:,:,9*local:9*local+6]*weight[:,None,None]).reshape(-1,6)
            A1[sl,6*index:6*index+6]=(W1[:,:,9*local:9*local+6]*weight[:,None,None]).reshape(-1,6)
            D1[sl,3*index:3*index+3]=(W1[:,:,9*local+6:9*local+9]*weight[:,None,None]).reshape(-1,3)
        points,_=common.element.bilinear(plane,xy[ids][[0,3,15,12]])
        for kind in targets:targets[kind][sl]=(common.target(kind,k,points[:,0])@candidate.CH.T*weight[:,None]).reshape(-1)
    M,area=common.h1tool.h1_matrix(xy,cells,coords,free,order)
    return dict(xy=xy,cells=cells,coords=coords,free=free,A0=A0,A1=A1,D1=D1,targets=targets,M=M,area=area,zero_c=zero_c,map_keys=keys)
