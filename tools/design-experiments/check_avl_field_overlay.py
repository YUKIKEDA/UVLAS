"""Predeclared basis and physical-overlay design experiment; no product code."""
from pathlib import Path
from decimal import Decimal,localcontext
import hashlib,json,math,time
import numpy as np
from scipy.linalg import block_diag,solve
import avl_field_overlay_exact as geometry
import check_vortex_core_load as segment
import check_avl_collinear_kernel as horseshoe
import check_avl_fixed_bands as old
import check_q9_plate_convergence as memory

ROOT=Path(__file__).resolve().parents[2];FOLDER=ROOT/'docs/verification'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def scalar(x):
    if isinstance(x,np.generic):return x.item()
    raise TypeError(type(x).__name__)
def write(p,obj):p.write_text(json.dumps(obj,indent=2,default=scalar,allow_nan=False)+'\n')
def array_digest(a):return hashlib.sha256(a.tobytes()).hexdigest()

def basis_trials():
    raw=json.loads((FOLDER/'avl-solve-sensitivity-runs.json').read_text());records=[]
    for run in raw['runs'][:2]:
        folder=ROOT/run['directory']/'out';m=run['m'];n=run['n'];size=m*n
        a=np.fromfile(folder/'us-matrix.bin',dtype='<f8').reshape((size,size),order='F')
        rows=np.loadtxt(folder/'uv-rows.txt',skiprows=2);g=np.loadtxt(folder/'us-vectors.txt',skiprows=1)[:,4]
        w=np.loadtxt(folder/'us-adjoints.txt',skiprows=1).reshape((22,size,5))[:,:,2].T
        before=[array_digest(x) for x in [a,rows,g,w]]
        p=block_diag(*([np.tril(np.ones((m,m)))]*n));d=block_diag(*([np.eye(m)-np.eye(m,k=-1)]*n))
        mu=p@g;ar=a@d;direct=np.empty_like(ar)
        points=rows[:,1:4];normal=rows[:,4:7]
        for j in range(size):
            fronta=rows[j,7:10];frontb=rows[j,10:13]
            if j%m<m-1:backa=rows[j+1,7:10];backb=rows[j+1,10:13]
            else:
                backa=fronta.copy();backb=frontb.copy()
                for v in [backa,backb]:v[0]=abs(v[1])*math.tan(math.radians(35))+1.5-.225*abs(v[1])
            vertices=np.array([fronta,frontb,backb,backa]);ends=np.roll(vertices,-1,axis=0)
            velocity=segment.kernel(points,vertices,ends,0.).sum(axis=1)*(4*np.pi*horseshoe.PI)
            if j%m==m-1:
                include=np.zeros((size,1),dtype=bool)
                bound,_=horseshoe.stable_bound(points,backa[None,:],backb[None,:],include)
                velocity+=(bound+horseshoe.stable_trail(points,backa[None,:],backb[None,:]))[:,0]
            direct[:,j]=np.sum(velocity*normal,axis=1)
        wr=d.T@w;q=w.T@g;qr=wr.T@mu;wrong=w.T@mu
        solved=solve(ar,-rows[:,6]);scale=max(np.max(abs(g)),1.)
        record=dict(case=run['case'],m=m,n=n,identity_error=float(max(np.max(abs(d@p-np.eye(size))),np.max(abs(p@d-np.eye(size))))),
            roundtrip_relative=float(np.max(abs(d@mu-g))/scale),kernel_error=float(np.max(abs(direct-ar))),
            kernel_eligible=bool(np.all(abs(direct-ar)<=1e-10+1e-10*abs(ar))),
            normalwash_difference=float(np.max(abs(a@g-ar@mu))),
            dual_output_difference=float(np.max(abs(q-qr))),negative_untransformed_weight=float(np.max(abs(q-wrong))),
            transformed_solve_relative=float(np.max(abs(solved-mu))/max(np.max(abs(mu)),1.)),
            condition1_horseshoe=float(np.linalg.cond(a,1)),condition1_ring=float(np.linalg.cond(ar,1)),
            inputs_unchanged=before==[array_digest(x) for x in [a,rows,g,w]])
        record['eligible']=all(record[k]<=1e-10 for k in ['identity_error','roundtrip_relative','normalwash_difference','dual_output_difference','transformed_solve_relative']) and record['kernel_eligible'] and record['negative_untransformed_weight']>1e-3 and record['inputs_unchanged']
        records.append(record)
    return records

def manufactured():
    D=geometry.D;domain=geometry.rectangle()
    def make(front,g=1.,spans=((0.,1.),)):
        return geometry.Field(domain,tuple(geometry.Strip(D(lo),D(hi),((D(front[0]),D(front[1])),),(g,)) for lo,hi in spans))
    constant=make((0,0));step=make((.25,0),2.);half=make((.5,0));diagonal=make((0,1))
    micro=make((.5+2**-40,0));split=make((.25,0),2.,((0,.3),(.3,1)))
    result,a,v=geometry.overlay([constant,step,split],observe=True)
    checks=dict(constants=abs(result['integrals'][0]-1)<=1e-10 and abs(result['integrals'][1]-1.5)<=1e-10,
        subdivision=np.max(abs(v[:,1]-v[:,2]))==0 and result['integral_preservation']<=1e-10)
    record,a,v=geometry.overlay([half,diagonal]);difference=geometry.metrics(a,v[:,0]-v[:,1]);inner=math.fsum(a*v[:,0]*v[:,1])
    checks['crossing_exact']=max(abs(difference['l1']-.25),abs(difference['l2']**2-.25),abs(inner-.375))<=1e-10
    micro_record,a,v=geometry.overlay([half,micro]);micro_values=geometry.metrics(a,v[:,0]-v[:,1])
    checks['micro_shift']=max(abs(micro_values['l1']-2**-40),abs(micro_values['l2']**2-2**-40),abs(micro_values['linf']-1))<=1e-15
    samples=(np.arange(16)+.5)/16;negative=float(np.sum((samples>=.5)!=(samples>=.5+2**-40)))
    checks['negative_midpoint']=negative==0 and micro_values['l1']>0
    triple,a,v=geometry.overlay([half,diagonal,constant]);d1=v[:,1]-v[:,0];d2=v[:,2]-v[:,1]
    triple_inner=math.fsum(a*d1*d2);checks['three_fields']=abs(triple_inner+.125)<=1e-10
    four,a,v=geometry.overlay([half,diagonal,constant,step]);checks['four_fields']=max(abs(x-y) for x,y in zip(four['integrals'],[.5,.5,1.,1.5]))<=1e-10 and abs(math.fsum(a*(v[:,1]-v[:,0])*(v[:,2]-v[:,1]))+.125)<=1e-10
    coefficient_sum,a,v=geometry.overlay([step,make((.25,0),.75),make((.25,0),1.25)])
    checks['coefficient_linearity']=np.max(abs(v[:,0]-v[:,1]-v[:,2]))==0
    # Identical rational straight line, with a non-dyadic span subdivision.
    boundary=D(float(1/3))
    exact_split=geometry.Field(domain,(geometry.Strip(D(0),boundary,((D(0),D(1)),),(1.,)),
        geometry.Strip(boundary,D(1),((D(0),D(1)),),(1.,))))
    same,a,v=geometry.overlay([diagonal,exact_split])
    same_difference=geometry.metrics(a,v[:,0]-v[:,1])
    checks['exact_collinearity']=same_difference['linf']==0
    # A shared binary endpoint with different slopes is preserved exactly.
    y=float(1/3);x=.4;left=geometry.line((.1,0.),(x,y));right=geometry.line((x,y),(.8,1.))
    shared=geometry.Field(domain,(geometry.Strip(D(0),D(y),(left,),(1.,)),geometry.Strip(D(y),D(1),(right,),(1.,))))
    shared_result,a,v=geometry.overlay([shared,shared])
    checks['shared_endpoint']=geometry.metrics(a,v[:,0]-v[:,1])['linf']==0 and shared_result['integral_preservation']<=1e-10
    obs,a,v=geometry.overlay([diagonal],observe=True);expected=[]
    for k in range(20):
        yl=k/20;yh=(k+1)/20
        for j in range(16):
            xl=j/16;xh=(j+1)/16
            expected.append(.5*(max(xh-yl,0)**2-max(xl-yl,0)**2-max(xh-yh,0)**2+max(xl-yh,0)**2))
    checks['fixed_cell_exact']=np.max(abs(np.array(obs['cell_integrals'][0])-expected))<=1e-10 and np.max(abs(np.array(obs['cell_areas'])-1/320))<=1e-10
    buffer=np.full(320,17.);before=array_digest(buffer);field_before=hashlib.sha256(repr(half).encode()).hexdigest();failures=[]
    invalid=[('nonfinite',lambda:make((.5,0),float('nan'))),
        ('gap',lambda:make((.5,0),spans=((0,.4),(.5,1)))),('overlap',lambda:make((.5,0),spans=((0,.6),(.5,1)))),
        ('degenerate',lambda:make((.5,0),spans=((0,0),(0,1)))),
        ('reversed_fronts',lambda:geometry.Field(domain,(geometry.Strip(D(0),D(1),((D(.8),D(0)),(D(.2),D(0))),(1.,1.)),))),
        ('crossed_internal_fronts',lambda:geometry.Field(domain,(geometry.Strip(D(0),D(1),((D(.2),D(.6)),(D(.8),D(-.6))),(1.,1.)),))),
        ('basis',lambda:geometry.Field(domain,half.strips,'ring')),
        ('domain',lambda:geometry.overlay([half,geometry.Field(geometry.Domain('other',domain.pieces),half.strips)]))]
    for name,action in invalid:
        try:action();refused=False
        except ValueError:refused=True
        failures.append(dict(kind=name,refused=refused,output_unchanged=array_digest(buffer)==before,
            input_unchanged=hashlib.sha256(repr(half).encode()).hexdigest()==field_before))
    checks['failures']=all(all(r[k] for k in ['refused','output_unchanged','input_unchanged']) for r in failures)
    return dict(checks=checks,crossing_norms=difference,crossing_inner=inner,micro=micro_values,negative_midpoint=negative,
        triple_difference_inner=triple_inner,fixed_cell_error=float(np.max(abs(np.array(obs['cell_integrals'][0])-expected))),failures=failures)

def load_field(run,family):
    p=ROOT/run['directory'];rows=np.loadtxt(p/'uv-rows.txt',skiprows=2);forces=np.loadtxt(p/'uv-forces.txt',skiprows=2)
    finite=forces[forces[:,2]==0];domain=geometry.wing_domain(run['case']);strips=[];error=0.
    for j in sorted(set(finite[:,1])):
        f=finite[finite[:,1]==j];ids=f[:,0].astype(int)-1;a=rows[ids,7:10];b=rows[ids,10:13]
        if np.max(abs(a[:,1]-a[0,1])) or np.max(abs(b[:,1]-b[0,1])):raise ValueError('nonshared_strip_endpoints')
        for points in [a,b]:
            expected_z=np.abs(points[:,1])*(0. if run['case']=='FW-V01' else math.tan(math.radians(5)))
            error=max(error,float(np.max(abs(points[:,2]-expected_z))))
        if error>1e-11:raise ValueError('off_physical_wing')
        strengths=tuple(float(x) for x in rows[ids,13]/math.sin(math.pi/180))
        fronts=tuple(geometry.line(aa,bb) for aa,bb in zip(a,b))
        strips.append(geometry.Strip(geometry.D(a[0,1]),geometry.D(b[0,1]),fronts,strengths))
    field=geometry.Field(domain,tuple(strips))
    g=rows[:,13]/math.sin(math.pi/180);a=rows[:,7:10];b=rows[:,10:13];mid=(a+b)/2;point=finite[:,3:6]
    f=2*g[:,None]*np.cross(np.array([1.,0.,0.]),b-a)
    line=math.fsum(np.cross(mid-old.ORIGIN,f)[:,1])/(old.S*old.C)
    couple=math.fsum(np.cross(point-mid,f)[:,1])/(old.S*old.C)
    lift=2*math.fsum(float(s.high-s.low)*math.fsum(s.strengths) for s in strips)/old.S
    return field,dict(family=family,case=run['case'],m=run['m'],n=run['n'],lift=lift,moment_line=line,moment_couple=couple,
        moment_point=line+couple,geometry_z_error=error,
        inputs={name:dict(path=(p/name).relative_to(ROOT).as_posix(),sha256=digest(p/name)) for name in ['uv-rows.txt','uv-forces.txt']})

def trace_overlay(fields):
    ys=sorted({y for f in fields for s in f.strips for y in [s.low,s.high]});areas=[];values=[]
    with localcontext() as ctx:
        ctx.prec=80
        for lo,hi in zip(ys[:-1],ys[1:]):
            areas.append(float(hi-lo));values.append([math.fsum(f.active((lo+hi)/2).strengths) for f in fields])
    return np.array(areas),np.array(values)

def cosine(a,x,y):
    nx=geometry.metrics(a,x)['l2'];ny=geometry.metrics(a,y)['l2']
    return None if nx==0 or ny==0 else math.fsum(a*x*y)/(nx*ny)

def main():
    begin=time.monotonic();out=FOLDER/'avl-field-overlay-result.json'
    if out.exists():raise FileExistsError(out)
    ledger=json.loads((FOLDER/'avl-solve-sensitivity-result.json').read_text())['circulation_ledger']
    basis=[];manufacture=None;records=[];series=[];stop=None;checks={};summaries=[]
    def budget():
        if time.monotonic()-begin>600:raise TimeoutError('overlay_time_budget')
        if memory.memory_bytes()>512*1024**2:raise MemoryError('overlay_memory_budget')
    try:
        basis=basis_trials();manufacture=manufactured()
        if not all(r['eligible'] for r in basis) or not all(manufacture['checks'].values()):raise RuntimeError('basis_or_manufacture_failed')
        fields={};observed={}
        for family,name in [('original','avl-fixed-band-grid-runs.json'),('aligned','avl-band-aligned-grid-runs.json')]:
            raw=json.loads((FOLDER/name).read_text())
            for run in raw['runs']:
                if run['alpha']!=1:continue
                field,rec=load_field(run,family);key=(family,run['case'],run['m'],run['n'])
                fields[key]=field;obs,_,_=geometry.overlay([field],observe=True);budget()
                obs['cell_means']=(np.array(obs['cell_integrals'][0])/obs['cell_areas']).tolist();rec['observations']=obs;observed[key]=obs
                previous=next(r for r in ledger if all(r[k]==rec[k] for k in ['family','case','m','n']))
                rec['load_ledger_difference']=max(abs(rec['lift']-previous['lift']),abs(rec['moment_point']-previous['moment_y']))
                records.append(rec);write(ROOT/'build/analysis/avl-field-overlay-progress.json',dict(completed_meshes=len(records),last=key))
            print(json.dumps(dict(stage='observations',family=family,meshes=len(records),seconds=time.monotonic()-begin)),flush=True)
        for family in ['original','aligned']:
            ms=[4,8,16,32] if family=='original' else [3,6,12,24]
            ns=[12,24,48,96] if family=='original' else [20,40,80,160]
            for case in ['FW-V01','FW-V02']:
                for direction,grid in [('chord',[(m,ns[-1]) for m in ms]),('span',[(ms[-1],n) for n in ns]),('diagonal',list(zip(ms,ns)))]:
                    keys=[(family,case,m,n) for m,n in grid];seq=[fields[k] for k in keys];obs=[observed[k] for k in keys]
                    pairs=[];triples=[]
                    for j in range(3):
                        summary,a,v=geometry.overlay(seq[j:j+2]);summaries.append(summary)
                        ta,tv=trace_overlay(seq[j:j+2]);cell_area=np.array(obs[j+1]['cell_areas']);cell=np.array(obs[j+1]['cell_means'])-obs[j]['cell_means']
                        pairs.append(dict(field=geometry.metrics(a,v[:,1]-v[:,0]),trace=geometry.metrics(ta,tv[:,1]-tv[:,0]),
                            cell=geometry.metrics(cell_area,cell),overlay=summary));budget()
                    for j in range(2):
                        summary,a,v=geometry.overlay(seq[j:j+3]);summaries.append(summary)
                        ta,tv=trace_overlay(seq[j:j+3]);cell_area=np.array(obs[j+2]['cell_areas'])
                        c0=np.array(obs[j]['cell_means']);c1=np.array(obs[j+1]['cell_means']);c2=np.array(obs[j+2]['cell_means'])
                        triples.append(dict(field=cosine(a,v[:,1]-v[:,0],v[:,2]-v[:,1]),trace=cosine(ta,tv[:,1]-tv[:,0],tv[:,2]-tv[:,1]),
                            cell=cosine(cell_area,c1-c0,c2-c1),overlay=summary));budget()
                    ta,tv=trace_overlay([seq[-1]])
                    scales=dict(field=obs[-1]['norms'][0]['l2'],trace=geometry.metrics(ta,tv[:,0])['l2'],
                        cell=geometry.metrics(np.array(obs[-1]['cell_areas']),np.array(obs[-1]['cell_means']))['l2'])
                    diagnostics={k:geometry.diagnostic([p[k]['l2'] for p in pairs],[t[k] for t in triples],scales[k]) for k in scales}
                    for k in diagnostics:diagnostics[k]['relative_last_difference']=pairs[-1][k]['l2']/scales[k] if scales[k]>0 else None
                    moment={}
                    rr=[next(r for r in records if (r['family'],r['case'],r['m'],r['n'])==k) for k in keys]
                    for k in ['lift','moment_line','moment_couple','moment_point']:
                        vals=[r[k] for r in rr];moment[k]=dict(values=vals,diagnostic=old.estimate(vals,max(abs(vals[-1]),.1)))
                    series.append(dict(family=family,case=case,direction=direction,meshes=grid,pairs=pairs,triples=triples,diagnostics=diagnostics,loads=moment))
                    print(json.dumps(dict(stage='series',family=family,case=case,direction=direction,status={k:v['status'] for k,v in diagnostics.items()},seconds=time.monotonic()-begin)),flush=True)
        all_overlay=[r['observations'] for r in records]+summaries
        checks=dict(basis=all(r['eligible'] for r in basis),manufactured=all(manufacture['checks'].values()),
            forty_meshes=len(records)==40 and len(series)==12,
            area_and_integral_preservation=all(abs(r['area']-old.S)<=1e-10 and r['integral_preservation']<=1e-10 for r in all_overlay),
            original_loads_preserved=all(r['load_ledger_difference']<=1e-10 for r in records),
            resource_budget=time.monotonic()-begin<600 and memory.memory_bytes()<512*1024**2)
    except Exception as exc:stop=repr(exc);raise
    finally:
        result=dict(protocol='docs/verification/avl-field-overlay-experiment.md',protocol_sha256=digest(FOLDER/'avl-field-overlay-experiment.md'),
            geometry_addendum_sha256=digest(FOLDER/'avl-field-overlay-exact-geometry-addendum.md'),geometry_arithmetic='exact rational input geometry; float64 measures and coefficients',
            source_sha256=digest(Path(__file__)),dependency_sha256={Path(m.__file__).name:digest(Path(m.__file__)) for m in [geometry,segment,horseshoe,old,memory]},
            parent_sha256={name:digest(FOLDER/name) for name in ['avl-solve-sensitivity-result.json','avl-solve-sensitivity-runs.json','avl-fixed-band-grid-runs.json','avl-band-aligned-grid-runs.json']},
            basis=basis,manufactured=manufacture,meshes=records,series=series,checks=checks,
            methods_criteria_pass=bool(checks) and all(checks.values()),stop_reason=stop,reference_convergence_accepted=False,
            seconds=time.monotonic()-begin,peak_process_working_bytes=memory.memory_bytes())
        write(out,result)
        if out.stat().st_size>64*1024**2:raise RuntimeError('overlay_disk_budget')
    print(json.dumps(dict(checks=checks,seconds=result['seconds'],peak=result['peak_process_working_bytes'],bytes=out.stat().st_size),indent=2))

if __name__=='__main__':main()
