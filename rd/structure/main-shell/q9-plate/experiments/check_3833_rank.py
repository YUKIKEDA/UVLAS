"""Limited design experiment; not a UVLAS solver implementation.
Protocol: rd/structure/main-shell/q9-plate/records/3833-algebraic-experiment.md.
"""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import json
import platform
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss

nodes=np.array([[-1,-1],[1,-1],[1,1],[-1,1],[0,-1],[1,0],[0,1],[-1,0]],dtype=float)
width=length=1.0
thickness=0.1
young=1.0
poisson=0.3
density=1.0

def shape(x,y):
    n=np.zeros(8); nx=np.zeros(8); ny=np.zeros(8)
    for i,(a,b) in enumerate(nodes[:4]):
        n[i]=.25*(1+a*x)*(1+b*y)*(a*x+b*y-1)
        nx[i]=a*.25*(1+b*y)*(2*a*x+b*y)
        ny[i]=b*.25*(1+a*x)*(a*x+2*b*y)
    n[4:]=[.5*(1-x*x)*(1-y),.5*(1+x)*(1-y*y),.5*(1-x*x)*(1+y),.5*(1-x)*(1-y*y)]
    nx[4:]=[-x*(1-y),.5*(1-y*y),-x*(1+y),-.5*(1-y*y)]
    ny[4:]=[-.5*(1-x*x),-y*(1+x),.5*(1-x*x),-y*(1-x)]
    if max(abs(n.sum()-1),abs(nx.sum()),abs(ny.sum()))>1e-12:
        raise RuntimeError('Shape identity failed')
    return n,nx,ny

for i,(x,y) in enumerate(nodes):
    if np.max(np.abs(shape(x,y)[0]-np.eye(8)[i]))>1e-12:
        raise RuntimeError('Nodal interpolation failed')

mu=young/(2*(1+poisson))
lam=young*poisson/((1+poisson)*(1-2*poisson))
elastic=np.diag([2*mu]*3+[mu]*3)
elastic[:3,:3]+=lam
sqrt_elastic=np.linalg.cholesky(elastic).T

def observations(order):
    mass_rows=[]; strain_rows=[]
    gx,wx=leggauss(order[0]); gy,wy=leggauss(order[1]); gz,wz=leggauss(order[2])
    for x,vx in zip(gx,wx):
        for y,vy in zip(gy,wy):
            n,nx,ny=shape(x,y)
            for zz,vz in zip(gz,wz):
                z=thickness*zz/2
                weight=length*width*thickness/8*vx*vy*vz
                scalar=np.column_stack((n,z*n,.5*z*z*n)).reshape(24)
                mass_rows.append(np.sqrt(density*weight)*np.kron(scalar.reshape(1,24),np.eye(3)))
                dx=np.column_stack((nx,z*nx,.5*z*z*nx)).reshape(24)*2/length
                dy=np.column_stack((ny,z*ny,.5*z*z*ny)).reshape(24)*2/width
                dz=np.column_stack((np.zeros(8),n,z*n)).reshape(24)
                bmat=np.zeros((6,72))
                for i,(a,b,c) in enumerate(zip(dx,dy,dz)):
                    bmat[:,3*i+0]=[a,0,0,0,c,b]
                    bmat[:,3*i+1]=[0,b,0,c,0,a]
                    bmat[:,3*i+2]=[0,0,c,b,a,0]
                strain_rows.append(np.sqrt(weight)*sqrt_elastic@bmat)
    return np.vstack(mass_rows),np.vstack(strain_rows)

def rank_info(a):
    s=np.linalg.svd(a,compute_uv=False)
    ranks={str(tol):int(np.count_nonzero(s>tol*s[0])) for tol in [1e-11,1e-10,1e-9]}
    rank=ranks['1e-10']
    return dict(ranks=ranks,stable=len(set(ranks.values()))==1,largest=float(s[0]),
                smallest_nonzero=float(s[rank-1]) if rank else None,
                largest_discarded=float(s[rank]) if rank<len(s) else 0.0)

rigid=np.zeros((72,6))
for i,(x,y) in enumerate(nodes):
    r=np.array([length*x/2,width*y/2,0.]); d=np.array([0.,0.,1.])
    rigid[9*i:9*i+3,:3]=np.eye(3)
    for j,omega in enumerate(np.eye(3)):
        rigid[9*i:9*i+3,3+j]=np.cross(omega,r)
        rigid[9*i+3:9*i+6,3+j]=np.cross(omega,d)

counter=np.zeros(72)
for i in range(8):
    counter[9*i:9*i+3]=-thickness**2/12*np.array([0.,0.,1.])
    counter[9*i+6:9*i+9]=2*np.array([0.,0.,1.])
exact_energy=density*length*width*thickness**5/360

results={}
for order in [(5,5,5),(3,3,3),(3,3,2),(2,2,2)]:
    m,b=observations(order)
    info=dict(mass=rank_info(m),stiffness=rank_info(b),
              rigid_residual=float(np.linalg.norm(b@rigid)/(np.linalg.norm(b)*np.linalg.norm(rigid))),
              counter_mass_energy=float(.5*np.linalg.norm(m@counter)**2))
    results['x'.join(map(str,order))]=info

checks={
    'rigid_basis_rank6':rank_info(rigid)['ranks']['1e-10']==6,
    'mass_full_rank72':results['5x5x5']['mass']['ranks']['1e-10']==72,
    'mass_3x3x3_rank72':results['3x3x3']['mass']['ranks']['1e-10']==72,
    'mass_3x3x2_rank48':results['3x3x2']['mass']['ranks']['1e-10']==48,
    'stiffness_full_rank66':results['5x5x5']['stiffness']['ranks']['1e-10']==66,
    'stiffness_2x2x2_rank_at_most48':results['2x2x2']['stiffness']['ranks']['1e-10']<=48,
    'all_rigid_residuals':all(v['rigid_residual']<=1e-12 for v in results.values()),
    'full_counter_energy':abs(results['5x5x5']['counter_mass_energy']/exact_energy-1)<=1e-12,
    'three_point_counter_energy':abs(results['3x3x3']['counter_mass_energy']/exact_energy-1)<=1e-12,
    'two_point_counter_zero':results['3x3x2']['counter_mass_energy']<=1e-20,
    'all_ranks_threshold_stable':all(v[k]['stable'] for v in results.values() for k in ['mass','stiffness']),
}
report=dict(protocol='rd/structure/main-shell/q9-plate/records/3833-algebraic-experiment.md',
            python=platform.python_version(),numpy=np.__version__,platform=platform.platform(),
            inputs=dict(a=length,b=width,h=thickness,Y=young,nu=poisson,rho=density),
            rank_method='SVD of weighted observation, relative cutoffs 1e-11 / 1e-10 / 1e-9',
            expected_counter_energy=exact_energy,checks=checks,results=results,
            all_protocol_checks_pass=all(checks.values()))
out=rd_paths.ROOT/'build'/'analysis'/'3833-rank-result.json'
out.parent.mkdir(parents=True,exist_ok=True)
out.write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
print(json.dumps(report,indent=2))
if not all(checks.values()):
    raise SystemExit(1)
