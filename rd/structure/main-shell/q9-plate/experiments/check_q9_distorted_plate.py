"""Prospective distorted planar mesh experiment; no production solver changes."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
from functools import lru_cache
import hashlib,itertools,json,math,platform,time
import numpy as np
import scipy
from numpy.polynomial.legendre import leggauss
from scipy.sparse import coo_matrix
import check_q9_independent_load as element
import check_q9_plate_convergence as plate
import check_q9_precision_diagnostic as precision
import check_q9_curved_load as geometry
ROOT=rd_paths.ROOT;F=rd_paths.RECORDS
PROTOCOL=F/'q9-distorted-plate-experiment.md'
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def mesh(n,amplitude=.2):
    u=np.linspace(0,1,n+1);vertices=np.array([[[x+(amplitude/n*(-1)**(i+j)*np.sin(np.pi*x)*np.sin(np.pi*y) if 0<i<n and 0<j<n else 0),y] for i,x in enumerate(u)] for j,y in enumerate(u)])
    nodes=np.full(((2*n+1)**2,2),np.nan);cells=plate.connectivity(n)
    for i,j,ids,_ in cells:
        corners=vertices[[j,j,j+1,j+1],[i,i+1,i+1,i]]
        for index,(x,y) in zip(ids,element.nodes):
            weights=np.array([(1-x)*(1-y),(1+x)*(1-y),(1+x)*(1+y),(1-x)*(1+y)])/4
            position=weights@corners
            if np.all(np.isfinite(nodes[index])):assert np.max(abs(nodes[index]-position))<1e-14
            nodes[index]=position
    return vertices,nodes,cells

def configuration(nodes):
    q=np.zeros((9,3,3));q[:,0,:2]=nodes;q[:,1,2]=1
    return q.reshape(81)

def determinant_bounds(nodes):
    determinants=[]
    for x,y in element.nodes[:4]:
        _,dx,dy=element.shape(x,y);j=nodes.T@np.column_stack((dx,dy));determinants.append(float(np.linalg.det(j)))
    return min(determinants),max(determinants)

@lru_cache(None)
def table(order,h):
    p,w=leggauss(order);z,wz=leggauss(3)
    points=[(x,y,t) for x,y,t in itertools.product(p,p,z)]
    weights=np.array([a*b*c for a,b,c in itertools.product(w,w,wz)])
    scalar,gradient=geometry.shape_tables(points,h)
    return scalar,gradient,weights

def stiffness(nodes,h,order):
    _,g,weights=table(order,h);q=configuration(nodes).reshape(27,3)
    jac=np.einsum('ia,pib->pab',q,g);det=np.linalg.det(jac)
    if np.min(det)<=0:raise ValueError('nonpositive volume Jacobian')
    grad=np.linalg.solve(jac.transpose(0,2,1),g.transpose(0,2,1)).transpose(0,2,1)
    b=np.zeros((len(g),6,81))
    b[:,0,0::3]=grad[:,:,0];b[:,1,1::3]=grad[:,:,1];b[:,2,2::3]=grad[:,:,2]
    b[:,3,1::3]=grad[:,:,2];b[:,3,2::3]=grad[:,:,1]
    b[:,4,0::3]=grad[:,:,2];b[:,4,2::3]=grad[:,:,0]
    b[:,5,0::3]=grad[:,:,1];b[:,5,1::3]=grad[:,:,0]
    rows=(np.sqrt(weights*det)[:,None,None]*np.einsum('ab,pbc->pac',element.chol,b)).reshape(-1,81)
    return rows.T@rows

@lru_cache(None)
def surface(order):
    p,w=leggauss(order);points=list(itertools.product(p,p))
    values=[element.shape(x,y) for x,y in points]
    return np.array([v[0] for v in values]),np.array([np.column_stack(v[1:]) for v in values]),np.array([a*b for a,b in itertools.product(w,w)])

def surface_data(nodes,order):
    shape,g,weights=surface(order);positions=shape@nodes
    jac=np.einsum('ia,pib->pab',nodes,g)
    return shape,positions,weights*np.linalg.det(jac)

def assemble(n,h,ref,nodes,cells,order,guard):
    values=[];f=np.zeros(len(nodes)*9)
    for _,_,ids,dofs in cells:
        guard();local=nodes[ids];values.append(stiffness(local,h,order).reshape(-1))
        shape,x,w=surface_data(local,order)
        force=shape.T@(ref['pressure']*np.sin(np.pi*x[:,0])*np.sin(np.pi*x[:,1])*w)
        v=np.zeros((9,9));v[:,2]=force;v[:,8]=h*h/8*force;np.add.at(f,dofs,v.reshape(81))
    rows=np.concatenate([np.repeat(d,81) for _,_,_,d in cells]);cols=np.concatenate([np.tile(d,81) for _,_,_,d in cells])
    k=coo_matrix((np.concatenate(values),(rows,cols)),shape=(len(f),len(f))).tocsr();k.eliminate_zeros()
    return k,f

def center_value(n,vertices,nodes,cells,q):
    j=n//2;t=n*.5-j;line=(1-t)*vertices[j,:,0]+t*vertices[j+1,:,0]
    i=min(n-1,max(0,int(np.searchsorted(line,.5,side='right')-1)))
    xi=2*(.5-line[i])/(line[i+1]-line[i])-1;eta=2*t-1
    _,_,ids,_=cells[j*n+i];shape,_,_=element.shape(xi,eta)
    assert np.max(abs(shape@nodes[ids]-.5))<=1e-12
    return float(shape@q[9*ids+2])

def observables(n,vertices,nodes,cells,q,f,ref):
    error=0.;den=0.
    for _,_,ids,_ in cells:
        shape,x,w=surface_data(nodes[ids],9)
        exact=ref['center']*np.sin(np.pi*x[:,0])*np.sin(np.pi*x[:,1])
        error+=float(w@((shape@q[9*ids+2]-exact)**2));den+=float(w@(exact**2))
    center=center_value(n,vertices,nodes,cells,q);compliance=math.fsum(float(x)*float(y) for x,y in zip(f,q))
    return dict(center=center,compliance=compliance,center_relative=abs(center/ref['center']-1),l2_relative=math.sqrt(error/den),compliance_relative=abs(compliance/ref['compliance']-1))

def patch(nodes,cells,n):
    error=0.;wrong=0.;area=0.;minimum=math.inf
    for _,_,ids,_ in cells:
        xy=nodes[ids];lo,_=determinant_bounds(xy);minimum=min(minimum,4*n*n*lo)
        assert lo>0 and 4*n*n*lo>=.2
        shape,g,weights=surface(5);j=np.einsum('ia,pib->pab',xy,g);grad=np.linalg.solve(j.transpose(0,2,1),g.transpose(0,2,1)).transpose(0,2,1)
        field=xy[:,0]**2;x=shape@xy;expected=np.column_stack((2*x[:,0],np.zeros(len(x))))
        error=max(error,float(np.max(abs(shape@field-x[:,0]**2))),float(np.max(abs(np.einsum('i,pia->pa',field,grad)-expected))))
        wrong=max(wrong,float(np.max(abs(np.einsum('i,pia->pa',field,g*2*n)-expected))))
        area+=float(weights@np.linalg.det(j))
    return dict(minimum_scaled_determinant=minimum,area=area,area_error=abs(area-1),polynomial_error=error,wrong_jacobian_error=wrong)

def run(n,h,out,stage_start):
    start=time.monotonic();deadline=min(start+240,stage_start+1200)
    def guard():
        if time.monotonic()>deadline:raise TimeoutError('case240_or_stage1200_seconds')
        if plate.memory_bytes()>2*1024**3:raise MemoryError('process_peak_2GiB')
        if sum(p.stat().st_size for p in out.rglob('*') if p.is_file())>64*1024**2:raise RuntimeError('output_64MiB')
    vertices,nodes,cells=mesh(n);geo=patch(nodes,cells,n)
    assert geo['polynomial_error']<=1e-11 and geo['area_error']<=1e-12 and geo['wrong_jacobian_error']>1e-3
    ref=plate.reference(h);free=precision.free_dofs(n);matrices={};loads={};solves={};solutions={};records={}
    for order in [7,9]:
        k,f=assemble(n,h,ref,nodes,cells,order,guard);matrices[order]=k;loads[order]=f
    common_scale=1/np.sqrt(matrices[7].diagonal());from scipy.sparse import diags
    d=diags(common_scale);base=d@matrices[7]@d;difference=d@(matrices[9]-matrices[7])@d
    quadrature=float(np.linalg.norm(difference.data)/np.linalg.norm(base.data));load_difference=float(np.linalg.norm(loads[9]-loads[7])/np.linalg.norm(loads[9]))
    symmetry=max(float(np.linalg.norm((k-k.T).data)/np.linalg.norm(k.data)) for k in matrices.values())
    assembly_pass=quadrature<=1e-8 and load_difference<=1e-8 and symmetry<=1e-12 and all(np.all(k.diagonal()>0) for k in matrices.values())
    if assembly_pass:
        for order in [7,9]:
            guard();a,b,scale,lu=precision.scaled_system(matrices[order],loads[order],free)
            solves[order]={}
            for digits in [70,90]:
                guard();r,q=precision.multiprecision(a,b,scale,lu,free,len(loads[order]),digits,deadline)
                r['observables']=observables(n,vertices,nodes,cells,q,loads[order],ref)
                solves[order][digits]=r;solutions[order,digits]=q
            records[str(order)]=precision.difference(solutions[order,70],solutions[order,90],solves[order][70]['observables'],solves[order][90]['observables'])
        integration_response=precision.difference(solutions[7,90],solutions[9,90],solves[7][90]['observables'],solves[9][90]['observables'])
        eligible=all(r['eligible'] for group in solves.values() for r in group.values()) and all(v<=1e-8 for group in records.values() for v in group.values()) and all(v<=1e-5 for v in integration_response.values())
    else:eligible=False;integration_response=None
    file=out/f'n{n}-h{h:g}.npz'
    np.savez_compressed(file,nodes=nodes,vertices=vertices,f7=loads[7],f9=loads[9],**{f'q{o}_{p}':q for (o,p),q in solutions.items()})
    guard()
    return dict(n=n,h=h,reference=ref,geometry=geo,stiffness_7_vs_9=quadrature,load_7_vs_9=load_difference,symmetry=symmetry,assembly_pass=bool(assembly_pass),solves=solves,
        precision_comparisons=records,quadrature_response=integration_response,eligible=bool(eligible),observables=solves[9][90]['observables'] if assembly_pass else None,
        seconds=time.monotonic()-start,peak_bytes=plate.memory_bytes(),artifact=file.relative_to(ROOT).as_posix(),artifact_sha256=digest(file))

def controls():
    vertices,nodes,cells=mesh(2,0);h=.01;ref=plate.reference(h)
    k,f=assemble(2,h,ref,nodes,cells,7,lambda:None)
    original=precision.assemble(plate.element_stiffness(2,h),cells,len(f));from scipy.sparse import diags
    s=diags(1/np.sqrt(original.diagonal()));a=s@original@s;delta=s@(k-original)@s
    matrix_error=float(np.linalg.norm(delta.data)/np.linalg.norm(a.data));force_error=float(np.linalg.norm(f-plate.load(2,h,ref['pressure'],cells,7))/np.linalg.norm(f))
    inverted=nodes[cells[0][2]].copy();inverted[:,0]*=-1;rejected=determinant_bounds(inverted)[0]<=0
    assert matrix_error<=1e-10 and force_error<=1e-10 and rejected
    return dict(rectangular_stiffness_difference=matrix_error,rectangular_load_difference=force_error,inverted_geometry_rejected=bool(rejected),all_criteria_pass=True)

def main():
    target=F/'q9-distorted-plate-result.json'
    if target.exists():raise FileExistsError(target)
    out=ROOT/'build/analysis'/('q9-distorted-plate-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir()
    start=time.monotonic();cases=[];stopped=None;control=None;checks=[]
    try:
        control=controls()
        for h in [.01,.001]:
            for n in [2,4,8,16]:
                result=run(n,h,out,start);cases.append(result)
                print(json.dumps(dict(h=h,n=n,eligible=result['eligible'],observables=result['observables'],seconds=result['seconds'])),flush=True)
        for h in [.01,.001]:
            a=next(r for r in cases if r['h']==h and r['n']==8);b=next(r for r in cases if r['h']==h and r['n']==16)
            rate=math.log2(a['observables']['l2_relative']/b['observables']['l2_relative']) if a['observables'] and b['observables'] else None
            checks.append(dict(h=h,accuracy_pass=bool(b['eligible'] and all(b['observables'][k]<=.01 for k in ['center_relative','l2_relative','compliance_relative'])),last_l2_order=rate,
                order_pass=bool(a['eligible'] and b['eligible'] and (b['observables']['l2_relative']<=1e-8 or rate>=1.5))))
        mid=next(r for r in cases if r['h']==.01 and r['n']==8);thin=next(r for r in cases if r['h']==.001 and r['n']==8)
        uniform=bool(mid['eligible'] and thin['eligible'] and thin['observables']['l2_relative']<=2*mid['observables']['l2_relative']+1e-8)
    except Exception as exc:stopped=repr(exc);uniform=False
    data=dict(protocol=PROTOCOL.relative_to(ROOT).as_posix(),protocol_sha256=digest(PROTOCOL),source_sha256=digest(__file__),
        dependency_sha256={Path(m.__file__).name:digest(m.__file__) for m in [element,plate,precision,geometry]},python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
        controls=control,cases=cases,thickness_checks=checks,thin_uniformity_pass=uniform,stop_reason=stopped,seconds=time.monotonic()-start,peak_bytes=plate.memory_bytes(),
        all_criteria_pass=bool(stopped is None and len(cases)==8 and all(r['eligible'] for r in cases) and all(r['accuracy_pass'] and r['order_pass'] for r in checks) and uniform),
        product_shell_adopted=False,limitation='Flat isotropic linear, this nonnested distortion family and quadratic thickness; not general curved-shell or exact 3D validation')
    target.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8');print(json.dumps(dict(all_criteria_pass=data['all_criteria_pass'],checks=checks,stop_reason=stopped)))
    if not data['all_criteria_pass']:raise SystemExit(1)
if __name__=='__main__':main()
