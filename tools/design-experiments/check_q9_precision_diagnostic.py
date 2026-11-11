"""Separate precision experiment, not a production solver.
Protocol: docs/verification/q9-precision-diagnostic-experiment.md, sections 1-7.
"""
import hashlib
import json
import math
import platform
import time
from decimal import Decimal, localcontext
from fractions import Fraction as F
from pathlib import Path
import numpy as np
import scipy
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import splu
from numpy.polynomial.legendre import leggauss
import check_q9_plate_convergence as plate
import check_q9_independent_load as element


def assemble(local, cells, size):
    rows=np.concatenate([np.repeat(d,81) for _,_,_,d in cells])
    cols=np.concatenate([np.tile(d,81) for _,_,_,d in cells])
    k=coo_matrix((np.tile(local.reshape(-1),len(cells)),(rows,cols)),shape=(size,size)).tocsr()
    k.eliminate_zeros()
    return k


def free_dofs(n):
    side=2*n+1;fixed=[]
    for j in range(side):
        for i in range(side):
            components=set()
            if i in (0,side-1):components.update([1,2])
            if j in (0,side-1):components.update([0,2])
            fixed.extend(9*(j*side+i)+3*t+c for t in range(3) for c in components)
    return np.setdiff1d(np.arange(9*side*side),fixed)


def scaled_system(k,f,free):
    kff=k[free,:][:,free]
    scale=1/np.sqrt(kff.diagonal())
    a=(diags(scale)@kff@diags(scale)).tocsr()
    return a,scale*f[free],scale,splu(a.tocsc())


def standard(a,b,scale,lu,free,size):
    y=lu.solve(b);history=[]
    for step in range(4):
        r=plate.careful_residual(a,y,b)
        relative=float(np.linalg.norm(r)/np.linalg.norm(b))
        history.append(dict(corrections=step,fsum_residual=relative,
            float64_residual=float(np.linalg.norm(b-a@y)/np.linalg.norm(b))))
        if relative<=1e-9:break
        if step<3:y+=lu.solve(r)
    q=np.zeros(size);q[free]=scale*y
    return dict(eligible=bool(relative<=1e-9 and np.all(np.isfinite(q))),history=history),q


def multiprecision(a,b,scale,lu,free,size,digits,deadline):
    # A,b are exact conversions of the fixed float64 inputs, not a new matrix.
    with localcontext() as ctx:
        ctx.prec=digits
        ad=[Decimal.from_float(float(v)) for v in a.data]
        bd=[Decimal.from_float(float(v)) for v in b]
        y=[Decimal.from_float(float(v)) for v in lu.solve(b)]
        initial=list(y);history=[]
        normb=float(np.linalg.norm(b))
        normainf=float(np.max(np.asarray(abs(a).sum(axis=1)).reshape(-1)))
        for step in range(6):
            if time.monotonic()>deadline:raise TimeoutError('120 second case budget')
            residual=[];component=0.
            for i,bi in enumerate(bd):
                lo,hi=a.indptr[i:i+2]
                products=[ad[k]*y[a.indices[k]] for k in range(lo,hi)]
                ri=bi-sum(products,Decimal(0))
                denominator=abs(bi)+sum((abs(v) for v in products),Decimal(0))
                ci=0. if denominator==0 and ri==0 else (math.inf if denominator==0 else float(abs(ri)/denominator))
                component=max(component,ci);residual.append(float(ri))
            residual=np.array(residual)
            relative=float(np.linalg.norm(residual)/normb)
            denom=normainf*max(abs(float(v)) for v in y)+float(np.max(np.abs(b)))
            backward=float(np.max(np.abs(residual))/denom)
            history.append(dict(corrections=step,relative_residual=relative,
                backward_normwise=backward,backward_componentwise=component))
            if relative<=1e-9 and backward<=1e-14 and component<=1e-14:break
            if step<5:
                delta=lu.solve(residual)
                y=[v+Decimal.from_float(float(d)) for v,d in zip(y,delta)]
        q=np.zeros(size);q0=np.zeros(size)
        for j,dof in enumerate(free):
            factor=Decimal.from_float(float(scale[j]))
            q[dof]=float(factor*y[j]);q0[dof]=float(factor*initial[j])
        projected=np.array([float(v) for v in y])
        projection_fsum=float(np.linalg.norm(plate.careful_residual(a,projected,b))/normb)
        projection64=float(np.linalg.norm(b-a@projected)/normb)
        groups={}
        for kind,name in enumerate(['position','first_thickness_derivative','second_thickness_derivative']):
            ids=free[(free%9)//3==kind]
            groups[name]=dict(max_abs_initial=float(np.max(np.abs(q0[ids]))),
                max_abs_correction=float(np.max(np.abs(q[ids]-q0[ids]))))
        serialized='\n'.join(str(v) for v in y).encode('ascii')
    record=dict(decimal_digits=digits,eligible=bool(relative<=1e-9 and backward<=1e-14 and component<=1e-14 and np.all(np.isfinite(q))),
        history=history,float64_projected_fsum_residual=projection_fsum,
        float64_projected_standard_residual=projection64,physical_correction_groups=groups,
        retained_scaled_solution_sha256=hashlib.sha256(serialized).hexdigest())
    return record,q


def poly_derivative(p):return [F(i)*p[i] for i in range(1,len(p))] or [F(0)]


def poly_integral_product(p,q,half):
    total=F(0)
    for i,a in enumerate(p):
        for j,b in enumerate(q):
            degree=i+j
            if degree%2==0:total+=a*b*2*half**(degree+1)/F(degree+1)
    return total


def rational_stiffness(n):
    # Tensor polynomial moments: independent of quadrature and element B.
    h=F(1,1000);nu=F(3,10)
    lam=nu/((1+nu)*(1-2*nu));mu=1/(2*(1+nu))
    polys={-1:[F(0),-F(1,2),F(1,2)],0:[F(1),F(0),-F(1)],1:[F(0),F(1,2),F(1,2)]}
    thickness=[[F(1)],[F(0),F(1)],[F(0),F(0),F(1,2)]]
    nodes=[(-1,-1),(1,-1),(1,1),(-1,1),(0,-1),(1,0),(0,1),(-1,0),(0,0)]
    scalars=[(x,y,t) for x,y in nodes for t in range(3)]
    grams={}
    for i in [-1,0,1]:
        for j in [-1,0,1]:
            for di in [0,1]:
                for dj in [0,1]:
                    grams[i,j,di,dj]=poly_integral_product(poly_derivative(polys[i]) if di else polys[i],poly_derivative(polys[j]) if dj else polys[j],F(1))
    zg={}
    for i in range(3):
        for j in range(3):
            for di in [0,1]:
                for dj in [0,1]:
                    zg[i,j,di,dj]=poly_integral_product(poly_derivative(thickness[i]) if di else thickness[i],poly_derivative(thickness[j]) if dj else thickness[j],h/2)
    result=np.zeros((81,81))
    for i,(x,y,t) in enumerate(scalars):
        for j,(xx,yy,tt) in enumerate(scalars):
            gradient={}
            for a in range(3):
                for b in range(3):
                    gradient[a,b]=F((2*n)**(int(a<2)+int(b<2)),4*n*n)*grams[x,xx,int(a==0),int(b==0)]*grams[y,yy,int(a==1),int(b==1)]*zg[t,tt,int(a==2),int(b==2)]
            trace=sum(gradient[c,c] for c in range(3))
            for a in range(3):
                for b in range(3):
                    result[3*i+a,3*j+b]=float(lam*gradient[a,b]+mu*gradient[b,a]+(mu*trace if a==b else F(0)))
    return result


def observables(n,q,f,ref,cells):
    side=2*n+1;center=float(q[9*(n*side+n)+2]);err=0.;den=0.
    points,weights=leggauss(7)
    for i,j,ids,_ in cells:
        w=q[9*ids+2]
        for x,wx in zip(points,weights):
            for y,wy in zip(points,weights):
                shape,_,_=element.shape(x,y)
                expected=ref['center']*np.sin(np.pi*(i+(x+1)/2)/n)*np.sin(np.pi*(j+(y+1)/2)/n)
                weight=wx*wy/(4*n*n)
                err+=(shape@w-expected)**2*weight;den+=expected**2*weight
    compliance=float(f@q)
    return dict(center=center,compliance=compliance,center_relative=float(abs(center/ref['center']-1)),
        l2_relative=float(np.sqrt(err/den)),compliance_relative=float(abs(compliance/ref['compliance']-1)))


def difference(q,p,obs,reference):
    result=dict(physical_solution_relative=float(np.max(np.abs(q-p))/np.max(np.abs(p))),
        center_relative=float(abs(obs['center']/reference['center']-1)),
        compliance_relative=float(abs(obs['compliance']/reference['compliance']-1)))
    for kind,name in enumerate(['position','first_thickness_derivative','second_thickness_derivative']):
        actual=q.reshape(-1,9)[:,3*kind:3*kind+3]
        expected=p.reshape(-1,9)[:,3*kind:3*kind+3]
        numerator=float(np.max(np.abs(actual-expected)));denominator=float(np.max(np.abs(expected)))
        result[name+'_relative']=numerator/denominator if denominator else (0. if numerator==0 else math.inf)
    return result


def run(n,h,ref):
    start=time.monotonic();deadline=start+120
    cells=plate.connectivity(n);size=9*(2*n+1)**2
    local=plate.element_stiffness(n,h);k=assemble(local,cells,size)
    f=plate.load(n,h,ref['pressure'],cells,7);f9=plate.load(n,h,ref['pressure'],cells,9)
    load_error=float(np.linalg.norm(f-f9)/np.linalg.norm(f9))
    symmetry=float(np.linalg.norm((k-k.T).data)/np.linalg.norm(k.data))
    record=dict(n=n,h=h,total_dofs=size,load_7_vs_9_relative=load_error,symmetry=symmetry)
    if load_error>1e-8 or symmetry>1e-12 or np.any(k.diagonal()<=0):
        return dict(**record,eligible=False,status='stopped_assembly')
    free=free_dofs(n);a,b,scale,lu=scaled_system(k,f,free)
    record['scaling_diagnostics']=dict(free_dofs=len(free),physical_diagonal_min=float(np.min(k.diagonal()[free])),
        physical_diagonal_max=float(np.max(k.diagonal()[free])),scale_min=float(np.min(scale)),scale_max=float(np.max(scale)),
        scaled_rhs_norm=float(np.linalg.norm(b)),scaled_rhs_max=float(np.max(np.abs(b))))
    normal,q=standard(a,b,scale,lu,free,size)
    normal['observables']=observables(n,q,f,ref,cells);record['float64_path']=normal
    eligible=normal['eligible'];selected=q;method='float64_fsum'
    if h==.001 and n in (8,16):
        records={};solutions={}
        for digits in [70,90]:
            precision,qp=multiprecision(a,b,scale,lu,free,size,digits,deadline)
            precision['observables']=observables(n,qp,f,ref,cells)
            records[str(digits)]=precision;solutions[digits]=qp
        comparison=difference(solutions[70],solutions[90],records['70']['observables'],records['90']['observables'])
        record['high_precision_paths']=records;record['70_vs_90']=comparison
        eligible=all(p['eligible'] for p in records.values()) and all(v<=1e-8 for v in comparison.values())
        exact=rational_stiffness(n)
        matrix_difference=element.matrix_error(exact,local)
        sensitivity=dict(element_matrix_difference=matrix_difference,eligible=False)
        if all(v<=1e-10 for v in matrix_difference.values()):
            kr=assemble(exact,cells,size);ar,br,sr,lur=scaled_system(kr,f,free)
            precision,qr=multiprecision(ar,br,sr,lur,free,size,90,deadline)
            ob=observables(n,qr,f,ref,cells)
            change=difference(qr,solutions[90],ob,records['90']['observables'])
            sensitivity.update(solve=precision,observables=ob,response_difference=change,
                eligible=bool(precision['eligible'] and all(v<=1e-6 for v in change.values())))
        record['rational_element_sensitivity']=sensitivity
        eligible=eligible and sensitivity['eligible']
        selected=solutions[90];method='decimal_90_retained_solution_float64_observables'
    record['selected_method']=method
    record['selected_observables']=observables(n,selected,f,ref,cells)
    record['seconds']=time.monotonic()-start
    record['peak_process_working_bytes']=plate.memory_bytes()
    record['eligible']=bool(eligible and record['seconds']<=120 and record['peak_process_working_bytes']<=4*1024**3)
    record['status']='evaluated' if record['eligible'] else 'ineligible_precision_or_budget'
    return record


def main():
    refs=[plate.reference(h) for h in [.1,.01,.001]];cases=[];stopped=None
    for ref in refs:
        for n in [2,4,8,16]:
            try:r=run(n,ref['h'],ref)
            except Exception as exc:r=dict(n=n,h=ref['h'],eligible=False,status='exception',error=str(exc))
            cases.append(r)
            print(json.dumps(dict(h=ref['h'],n=n,status=r['status'],seconds=r.get('seconds'))),flush=True)
            if r['status']=='exception' or r.get('seconds',0)>120 or r.get('peak_process_working_bytes',0)>4*1024**3:
                stopped=r;break
        if stopped:break
    checks=[]
    for ref in refs:
        fine=next((r for r in cases if r['h']==ref['h'] and r['n']==16),None)
        prev=next((r for r in cases if r['h']==ref['h'] and r['n']==8),None)
        fine_obs=fine.get('selected_observables') if fine else None
        prev_obs=prev.get('selected_observables') if prev else None
        rate=float(np.log2(prev_obs['l2_relative']/fine_obs['l2_relative'])) if fine_obs and prev_obs else None
        accuracy=bool(fine and fine['eligible'] and all(fine_obs[v]<=.01 for v in ['center_relative','l2_relative','compliance_relative']))
        order=bool(fine and prev and fine['eligible'] and prev['eligible'] and (fine_obs['l2_relative']<=1e-8 or rate>=1.5))
        checks.append(dict(h=ref['h'],finest_accuracy_pass=accuracy,last_l2_order=rate,order_pass=order))
    mid=next((c for c in cases if c['h']==.01 and c['n']==8),None)
    thin=next((c for c in cases if c['h']==.001 and c['n']==8),None)
    uniform=bool(mid and thin and mid['eligible'] and thin['eligible'] and thin['selected_observables']['l2_relative']<=2*mid['selected_observables']['l2_relative']+1e-8)
    modal=all(r['residual']<=1e-9 for r in refs) and refs[-1]['kl_relative']<=1e-3 and refs[-1]['kl_relative']<refs[-2]['kl_relative']
    passed=len(cases)==12 and all(c['eligible'] for c in cases) and all(c['finest_accuracy_pass'] and c['order_pass'] for c in checks) and uniform and modal
    report=dict(protocol='docs/verification/q9-precision-diagnostic-experiment.md',python=platform.python_version(),
        numpy=np.__version__,scipy=scipy.__version__,references=refs,cases=cases,thickness_checks=checks,
        thin_uniformity_pass=uniform,modal_limit_pass=bool(modal),all_criteria_pass=bool(passed),stopped=stopped,
        limitation='Fixed float64 assembled systems with retained Decimal solutions; rational element sensitivity is not exact global assembly or full 3D validation')
    path=Path(__file__).resolve().parents[2]/'build/analysis/q9-precision-diagnostic-result.json'
    path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print(json.dumps(dict(all_criteria_pass=passed,thickness_checks=checks)),flush=True)
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
