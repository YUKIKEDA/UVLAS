"""Limited linear plate experiment; not product code.
Protocol: docs/verification/q9-plate-convergence-experiment.md (sections 1-7).
"""
import ctypes
import json
import math
import platform
import time
from pathlib import Path
import numpy as np
import scipy
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import splu
from numpy.polynomial.legendre import leggauss
import check_q9_independent_load as element


def memory_bytes():
    if platform.system() != 'Windows':
        raise RuntimeError('Memory measurement currently specified for Windows only')
    class Counters(ctypes.Structure):
        _fields_ = [('cb', ctypes.c_ulong), ('faults', ctypes.c_ulong)] + [
            (n, ctypes.c_size_t) for n in ['peak_working', 'working', 'peak_paged',
            'paged', 'peak_nonpaged', 'nonpaged', 'pagefile', 'peak_pagefile']]
    value = Counters(); value.cb = ctypes.sizeof(value)
    handle = ctypes.windll.kernel32.GetCurrentProcess
    handle.restype = ctypes.c_void_p
    query = ctypes.windll.psapi.GetProcessMemoryInfo
    query.argtypes = [ctypes.c_void_p, ctypes.POINTER(Counters), ctypes.c_ulong]
    if not query(handle(), ctypes.byref(value), value.cb):
        raise RuntimeError('Cannot measure process memory')
    return int(value.peak_working)


def reference(h):
    # Independent modal energy: no element shape or stiffness helpers.
    young = 1.; poisson = .3; alpha = beta = np.pi
    mu = young/(2*(1+poisson))
    lam = young*poisson/((1+poisson)*(1-2*poisson))
    material = np.diag([2*mu]*3+[mu]*3)
    material[:3, :3] += lam
    k = np.zeros((4, 4))
    for t, weight in zip(*leggauss(5)):
        z = h*t/2
        g = np.array([[-alpha*z,0,0,0], [0,-beta*z,0,0], [0,0,0,z],
            [0,1,beta,.5*beta*z*z], [1,0,alpha,.5*alpha*z*z],
            [beta*z,alpha*z,0,0]])
        k += (g.T@material@g)*weight*h/8
    wkl = .001*h
    pressure = young*h**3/(12*(1-poisson**2))*(alpha**2+beta**2)**2*wkl
    f = pressure/4*np.array([0,0,1,h*h/8])
    scale = 1/np.sqrt(np.diag(k))
    a = scale[:, None]*k*scale[None, :]; rhs = scale*f
    s = scale*np.linalg.solve(a, rhs)
    residual = np.linalg.norm(a@(s/scale)-rhs)/np.linalg.norm(rhs)
    return dict(h=h,coefficients=s.tolist(),pressure=pressure,center=float(s[2]),
        compliance=float(f@s),kl_center=wkl,kl_relative=float(abs(s[2]/wkl-1)),
        residual=float(residual),condition=float(np.linalg.cond(a)))


def element_stiffness(n, h):
    k = np.zeros((81, 81))
    for x, wx in zip(*leggauss(5)):
        for y, wy in zip(*leggauss(5)):
            for t, wt in zip(*leggauss(5)):
                _, b = element.kinematics(x, y, h*t/2)
                b[[0,1,5], :] *= n
                b[3,2::3] *= n; b[4,2::3] *= n
                k += b.T@element.elastic@b*(wx*wy*wt*h/(8*n*n))
    return k


def connectivity(n):
    offsets = np.rint(element.nodes+1).astype(int)
    side = 2*n+1
    cells = []
    for j in range(n):
        for i in range(n):
            ids = (2*j+offsets[:,1])*side+2*i+offsets[:,0]
            cells.append((i,j,ids,(9*ids[:,None]+np.arange(9)).reshape(81)))
    return cells


def load(n, h, pressure, cells, order):
    f = np.zeros(9*(2*n+1)**2)
    points, weights = leggauss(order)
    for i,j,ids,dofs in cells:
        local = np.zeros((9,9))
        for x, wx in zip(points, weights):
            for y, wy in zip(points, weights):
                shape, _, _ = element.shape(x,y)
                px = (i+(x+1)/2)/n; py = (j+(y+1)/2)/n
                r = shape*pressure*np.sin(np.pi*px)*np.sin(np.pi*py)*wx*wy/(4*n*n)
                local[:,2] += r
                local[:,8] += h*h/8*r
        np.add.at(f,dofs,local.reshape(81))
    return f


def careful_residual(a, y, rhs):
    return np.array([b-math.fsum(float(v)*float(y[j])
        for v,j in zip(a.data[a.indptr[i]:a.indptr[i+1]],
                       a.indices[a.indptr[i]:a.indptr[i+1]]))
        for i,b in enumerate(rhs)])


def run(n,h,ref):
    start = time.monotonic(); cells = connectivity(n)
    local = element_stiffness(n,h); size = 9*(2*n+1)**2
    rows = np.concatenate([np.repeat(d,81) for _,_,_,d in cells])
    cols = np.concatenate([np.tile(d,81) for _,_,_,d in cells])
    values = np.tile(local.reshape(-1),len(cells))
    k = coo_matrix((values,(rows,cols)),shape=(size,size)).tocsr()
    del rows,cols,values
    norm = np.linalg.norm(k.data)
    symmetry = float(np.linalg.norm((k-k.T).data)/norm)
    f = load(n,h,ref['pressure'],cells,5)
    f7 = load(n,h,ref['pressure'],cells,7)
    load_error = float(np.linalg.norm(f-f7)/np.linalg.norm(f7))
    side = 2*n+1; fixed = []
    for j in range(side):
        for i in range(side):
            components = set()
            if i in (0,side-1): components.update([1,2])
            if j in (0,side-1): components.update([0,2])
            fixed.extend(9*(j*side+i)+3*t+c for t in range(3) for c in components)
    free = np.setdiff1d(np.arange(size),fixed)
    record = dict(h=h,n=n,total_dofs=size,free_dofs=len(free),symmetry=symmetry,
                  load_quadrature_relative=load_error)
    if symmetry>1e-12 or load_error>1e-8 or np.any(k.diagonal()<=0):
        record.update(status='stopped_assembly',criteria_pass=False)
        return record
    kff = k[free,:][:,free]
    scale = 1/np.sqrt(kff.diagonal())
    a = (diags(scale)@kff@diags(scale)).tocsr()
    rhs = scale*f[free]
    lu = splu(a.tocsc()); y = lu.solve(rhs); history=[]
    for step in range(4):
        residual = careful_residual(a,y,rhs)
        r = float(np.linalg.norm(residual)/np.linalg.norm(rhs))
        r64 = float(np.linalg.norm(rhs-a@y)/np.linalg.norm(rhs))
        history.append(dict(corrections=step,fsum_residual=r,float64_residual=r64))
        if r<=1e-9: break
        if step<3: y+=lu.solve(residual)
    q = np.zeros(size); q[free] = scale*y
    center = float(q[9*(n*side+n)+2])
    error_integral = 0.; reference_integral = 0.
    points,weights = leggauss(7)
    for i,j,ids,dofs in cells:
        w = q[9*ids+2]
        for x,wx in zip(points,weights):
            for yy,wy in zip(points,weights):
                shape,_,_ = element.shape(x,yy)
                exact = ref['center']*np.sin(np.pi*(i+(x+1)/2)/n)*np.sin(np.pi*(j+(yy+1)/2)/n)
                weight = wx*wy/(4*n*n)
                error_integral += (shape@w-exact)**2*weight
                reference_integral += exact**2*weight
    elapsed = time.monotonic()-start; memory = memory_bytes()
    eligible = r<=1e-9 and np.all(np.isfinite(q)) and elapsed<=120 and memory<=4*1024**3
    record.update(status='evaluated' if eligible else 'ineligible_solve_or_budget',
        eligible=bool(eligible),refinement_history=history,center=center,
        center_relative=float(abs(center/ref['center']-1)),
        l2_relative=float(np.sqrt(error_integral/reference_integral)),
        compliance_relative=float(abs(float(f@q)/ref['compliance']-1)),
        seconds=elapsed,peak_process_working_bytes=memory)
    return record


def main():
    refs = [reference(h) for h in [.1,.01,.001]]
    cases=[]; stopped=None
    for ref in refs:
        for n in [2,4,8,16]:
            try: result=run(n,ref['h'],ref)
            except Exception as exc:
                result=dict(h=ref['h'],n=n,status='exception',error=str(exc),eligible=False)
            cases.append(result)
            if result.get('seconds',0)>120 or result.get('peak_process_working_bytes',0)>4*1024**3 or result['status']=='exception':
                stopped=result;break
        if stopped:break
    checks=[]
    for ref in refs:
        series=[c for c in cases if c['h']==ref['h']]
        finest=next((c for c in series if c['n']==16),None)
        previous=next((c for c in series if c['n']==8),None)
        accuracy=bool(finest and finest.get('eligible') and all(finest[x]<=.01 for x in ['center_relative','l2_relative','compliance_relative']))
        rate=None if not finest or not previous or 'l2_relative' not in finest or 'l2_relative' not in previous else float(np.log2(previous['l2_relative']/finest['l2_relative']))
        order_ok=bool(finest and finest.get('eligible') and previous and previous.get('eligible') and (finest['l2_relative']<=1e-8 or (rate is not None and rate>=1.5)))
        checks.append(dict(h=ref['h'],finest_accuracy_pass=accuracy,last_l2_order=rate,order_pass=order_ok))
    mid=next((c for c in cases if c['h']==.01 and c['n']==8),None)
    thin=next((c for c in cases if c['h']==.001 and c['n']==8),None)
    uniform=bool(mid and thin and mid.get('eligible') and thin.get('eligible') and thin['l2_relative']<=2*mid['l2_relative']+1e-8)
    modal=bool(all(r['residual']<=1e-9 for r in refs) and refs[-1]['kl_relative']<=1e-3 and refs[-1]['kl_relative']<refs[-2]['kl_relative'])
    passed=all(c.get('eligible',False) for c in cases) and len(cases)==12 and all(c['finest_accuracy_pass'] and c['order_pass'] for c in checks) and uniform and modal
    report=dict(protocol='docs/verification/q9-plate-convergence-experiment.md',
        python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
        references=refs,cases=cases,thickness_checks=checks,
        thin_uniformity_pass=uniform,modal_kl_limit_pass=modal,
        all_criteria_pass=bool(passed),stopped=stopped,
        limitation='Flat isotropic linear rectangular mesh; modal reference restricts thickness to quadratic; not exact 3D elasticity')
    out=Path(__file__).resolve().parents[2]/'build/analysis/q9-plate-convergence-result.json'
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print(json.dumps(report,indent=2))
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
