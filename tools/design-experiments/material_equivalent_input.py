"""Explicit 3D elastic/effective input for design verification; SI/Mandel."""
from dataclasses import dataclass
from fractions import Fraction
import hashlib,json
import numpy as np
MB=np.zeros((6,3,3));MB[0,0,0]=MB[1,1,1]=MB[2,2,2]=1.
for k,(i,j) in enumerate([(1,2),(0,2),(0,1)],3):MB[k,i,j]=MB[k,j,i]=1/np.sqrt(2.)
TAGS=dict(units='SI',components='mandel_xx_yy_zz_yz_xz_xy',measure='green_second_pk')

def readonly(a):
    a=np.array(a,dtype=float,copy=True);a.flags.writeable=False;return a
def tensor(e):return np.einsum('...k,kij->...ij',e,MB)
def transform(A):return np.einsum('kij,...ia,lab,...jb->...kl',MB,A,MB,A,optimize=True)
def finite(a,shape,name):
    a=np.asarray(a,dtype=float)
    if a.shape!=shape or not np.isfinite(a).all():raise ValueError(name)
    return a
def positive(C):
    C=finite(C,(6,6),'stiffness_shape_or_finite')
    asym=float(np.linalg.norm(C-C.T)/max(np.linalg.norm(C),1e-300))
    if asym>1e-12:raise ValueError('nonsymmetric_stiffness')
    if np.min(np.diag(C))<=0:raise ValueError('nonpositive_stiffness_diagonal')
    C=(C+C.T)/2;scale=1/np.sqrt(np.diag(C));minimum=float(np.linalg.eigvalsh(C*scale[:,None]*scale[None,:]).min())
    if minimum<=1e-12:raise ValueError('nonpositive_or_unresolved_stiffness')
    return readonly(C),dict(symmetry_correction=asym/2,scaled_minimum=minimum)

@dataclass(frozen=True)
class Material:
    C:np.ndarray
    e0:np.ndarray
    rho:float
    diagnostics:dict

def prepare_material(request):
    if any(request.get(k)!=v for k,v in TAGS.items() if k!='components') or request.get('components') not in [TAGS['components'],'engineering_xx_yy_zz_yz_xz_xy']:raise ValueError('units_components_or_measure')
    kind=request.get('kind')
    if 'rho' not in request or 'e0' not in request:raise ValueError('missing_density_or_initial_strain')
    if kind=='isotropic':
        if 'E' not in request or 'nu' not in request:raise ValueError('missing_isotropic_constants')
        E=float(request['E']);nu=float(request['nu'])
        if not np.isfinite([E,nu]).all() or E<=0 or not -1<nu<.5:raise ValueError('invalid_isotropic_constants')
        mu=E/(2*(1+nu));lam=E*nu/((1+nu)*(1-2*nu));C=np.eye(6)*2*mu;C[:3,:3]+=lam
    elif kind=='orthotropic':
        keys=['E1','E2','E3','nu12','nu13','nu23','G23','G13','G12']
        if any(k not in request for k in keys):raise ValueError('missing_3d_constants')
        E1,E2,E3,n12,n13,n23,G23,G13,G12=[float(request[k]) for k in keys]
        if not np.isfinite([E1,E2,E3,n12,n13,n23,G23,G13,G12]).all() or min(E1,E2,E3,G23,G13,G12)<=0:raise ValueError('invalid_orthotropic_constants')
        H=np.diag([1/E1,1/E2,1/E3,1/G23,1/G13,1/G12]);H[0,1]=H[1,0]=-n12/E1;H[0,2]=H[2,0]=-n13/E1;H[1,2]=H[2,1]=-n23/E2
        scale=1/np.sqrt(np.diag(H))
        if np.linalg.eigvalsh(H*scale[:,None]*scale[None,:]).min()<=1e-12:raise ValueError('nonpositive_compliance')
        D=np.diag([1,1,1,np.sqrt(2),np.sqrt(2),np.sqrt(2)]);C=D@np.linalg.solve(H,np.eye(6))@D
    elif kind=='effective_3d':
        if 'C' not in request:raise ValueError('missing_3d_stiffness')
        C=finite(request['C'],(6,6),'stiffness_shape_or_finite')
        if request['components'].startswith('engineering'):
            D=np.diag([1,1,1,np.sqrt(2),np.sqrt(2),np.sqrt(2)]);C=D@C@D
    elif kind in ['abd','abd_shear']:raise ValueError('insufficient_3d_closure')
    else:raise ValueError('unsupported_material_kind')
    C,diagnostics=positive(C);rho=float(request['rho'])
    if not np.isfinite(rho) or rho<=0:raise ValueError('nonpositive_density')
    e0=finite(request['e0'],(6,),'initial_strain_shape_or_finite')
    if request['components'].startswith('engineering'):e0=e0/np.array([1,1,1,np.sqrt(2),np.sqrt(2),np.sqrt(2)])
    natural=float(np.linalg.eigvalsh(np.eye(3)+2*tensor(e0)).min())
    if natural<=1e-12:raise ValueError('invalid_natural_metric')
    return Material(C,readonly(e0),rho,dict(diagnostics,natural_metric_minimum=natural))

@dataclass(frozen=True)
class Layer:
    id:str
    lower:Fraction
    upper:Fraction
    material:Material
    fiber:np.ndarray

def prepare_stack(requests):
    if not requests:raise ValueError('empty_stack')
    out=[];cursor=Fraction(-1);ids=set()
    for request in requests:
        lo=Fraction(str(request['lower']));hi=Fraction(str(request['upper']));ident=request['id']
        if not isinstance(ident,str) or not ident or ident in ids:raise ValueError('invalid_or_duplicate_layer_id')
        if lo!=cursor or hi<=lo or hi>1:raise ValueError('layer_gap_overlap_or_order')
        fiber=finite(request['fiber_global'],(3,),'fiber_shape_or_finite')
        if np.linalg.norm(fiber)==0:raise ValueError('zero_fiber')
        out.append(Layer(ident,lo,hi,prepare_material(request['material']),readonly(fiber)));ids.add(ident);cursor=hi
    if cursor!=1:raise ValueError('incomplete_stack')
    return tuple(out)

def layer_at(stack,zeta,side=None):
    z=Fraction(str(zeta))
    if z< -1 or z>1:raise ValueError('outside_thickness')
    if side not in [None,'lower','upper']:raise ValueError('invalid_side')
    boundary=any(z==l.upper for l in stack[:-1])
    if side is not None and not boundary:raise ValueError('side_without_interface')
    if boundary and side=='lower':return next(l for l in stack if l.upper==z)
    for i,l in enumerate(stack):
        if l.lower<=z<l.upper or i==len(stack)-1 and z==l.upper:return l
    raise ValueError('no_layer')

def axes(J,fiber):
    J=np.asarray(J);normal=np.cross(J[:,:,0],J[:,:,1]);norm=np.linalg.norm(normal,axis=1)
    if np.any(norm<=0) or np.any(np.linalg.det(J)<=0):raise ValueError('invalid_reference_axes')
    normal=normal/norm[:,None];v=fiber-np.sum(fiber*normal,axis=1)[:,None]*normal;vn=np.linalg.norm(v,axis=1)
    if np.any(vn<=1e-12*np.linalg.norm(fiber)):raise ValueError('singular_fiber_projection')
    a=v/vn[:,None];return np.stack((a,np.cross(normal,a),normal),axis=2)

def point_fields(layer,J):
    A=axes(J,layer.fiber);T=transform(A);m=layer.material
    return T@m.C@T.transpose(0,2,1),np.einsum('pij,j->pi',T,m.e0),A,T

def identity(stack):
    data=[dict(id=l.id,lower=str(l.lower),upper=str(l.upper),C=l.material.C.tolist(),rho=l.material.rho,e0=l.material.e0.tolist(),fiber=l.fiber.tolist()) for l in stack]
    return hashlib.sha256(json.dumps(dict(tags=TAGS,layers=data),sort_keys=True,allow_nan=False,separators=(',',':')).encode()).hexdigest()

def section_moments(stack,h):
    if not np.isfinite(h) or h<=0:raise ValueError('nonpositive_thickness')
    C=[np.zeros((6,6)) for _ in range(5)];g=[np.zeros(6) for _ in range(3)];rho=np.zeros(5);constant=0.
    for layer in stack:
        lo=float(layer.lower)*h/2;hi=float(layer.upper)*h/2;m=layer.material
        A=axes(np.eye(3)[None],layer.fiber)[0];T=transform(A);Cg=T@m.C@T.T;e=T@m.e0
        moments=np.array([(hi**(k+1)-lo**(k+1))/(k+1) for k in range(5)])
        for k in range(5):C[k]+=moments[k]*Cg;rho[k]+=moments[k]*m.rho
        for k in range(3):g[k]+=moments[k]*Cg@e
        constant+=moments[0]*(e@Cg@e)/2
    return np.array(C),np.array(g),rho,constant

def verify_certificate(stack,h,certificate):
    if certificate.get('reference')!='mid_surface' or certificate.get('units')!='SI' or certificate.get('h')!=h:raise ValueError('certificate_reference_or_units')
    C,g,rho,w=section_moments(stack,h)
    for key,actual in [('C_moments',C),('initial_moments',g),('mass_moments',rho),('initial_constant',np.array(w))]:
        supplied=np.asarray(certificate.get(key),float)
        if supplied.shape!=actual.shape or not np.isfinite(supplied).all() or np.linalg.norm(supplied-actual)>1e-10*max(np.linalg.norm(actual),1e-300):raise ValueError('certificate_mismatch_'+key)
    return True
