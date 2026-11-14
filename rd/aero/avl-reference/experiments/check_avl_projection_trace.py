"""Predeclared manufactured projection and signed trace diagnostics only."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
from fractions import Fraction
import hashlib,json,math,time
import numpy as np
import avl_field_overlay_exact as geo
import check_avl_field_overlay as prior
import check_q9_plate_convergence as memory

ROOT=rd_paths.ROOT;FOLDER=rd_paths.RECORDS
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,obj):p.write_text(json.dumps(obj,indent=2,default=prior.scalar,allow_nan=False)+'\n')
def add(*polys):
    r={}
    for p in polys:
        for k,v in p.items():r[k]=r.get(k,Fraction(0))+v
    return {k:v for k,v in r.items() if v}
def scale(p,c):return {k:v*c for k,v in p.items()}
def multiply(p,q):
    r={}
    for (i,j),a in p.items():
        for (k,l),b in q.items():r[(i+k,j+l)]=r.get((i+k,j+l),Fraction(0))+a*b
    return {k:v for k,v in r.items() if v}
def target(piece,kind,cref):
    le,te=piece[2:];q={(1,0):Fraction(1),(0,0):-le[0],(0,1):-le[1]}
    c={(0,0):te[0]-le[0],(0,1):te[1]-le[1]}
    return q if kind==1 else scale(multiply(q,add(scale(c,2),scale(q,-1))),1/cref)

def integral(poly,left,right,lo,hi):
    """Exact monomial moments of a trapezoid, using rational coordinates."""
    total=Fraction(0)
    for (i,j),coefficient in poly.items():
        for edge,sign in [(right,1),(left,-1)]:
            for k in range(i+2):
                power=j+k+1
                term=Fraction(sign*math.comb(i+1,k),i+1)*edge[0]**(i+1-k)*edge[1]**k
                total+=coefficient*term*(hi**power-lo**power)/power
    return total

def gauss(poly,left,right,lo,hi):
    q,w=np.polynomial.legendre.leggauss(5);ys=float(lo)+(q+1)*float(hi-lo)/2
    terms=[]
    for y,wy in zip(ys,w):
        xl=float(left[0])+float(left[1])*y;xr=float(right[0])+float(right[1])*y
        xs=xl+(q+1)*(xr-xl)/2
        values=sum(float(c)*xs**i*y**j for (i,j),c in poly.items())
        terms.append(float(wy)*float(hi-lo)/2*(xr-xl)/2*float(np.dot(w,values)))
    return math.fsum(terms)

def project(field,kind,cref):
    strips=[];regions=[];orthogonal=0.;reference=[];reference_sq=[];means=[]
    for piece in field.domain.pieces:
        poly=target(piece,kind,cref)
        reference.append(float(integral(poly,piece[2],piece[3],piece[0],piece[1])))
        reference_sq.append(float(integral(multiply(poly,poly),piece[2],piece[3],piece[0],piece[1])))
    for s in field.strips:
        piece=field.domain.piece((s.low+s.high)/2);poly=target(piece,kind,cref);square=multiply(poly,poly)
        bounds=[piece[2]]+list(s.fronts)+[piece[3]];constant=[];local=[]
        for j,(left,right) in enumerate(zip(bounds[:-1],bounds[1:])):
            area=integral({(0,0):Fraction(1)},left,right,s.low,s.high)
            if area<0:raise ArithmeticError('negative_projection_area')
            i1=integral(poly,left,right,s.low,s.high);i2=integral(square,left,right,s.low,s.high)
            average=Fraction(0) if j==0 else i1/area
            if j>0:constant.append(float(average))
            local.append((j,left,right,area,i1,i2,poly))
        g=[constant[0]]+[constant[j]-constant[j-1] for j in range(1,len(constant))]
        strips.append(geo.Strip(s.low,s.high,s.fronts,tuple(g)))
        cumulative=[0.]+[math.fsum(g[:j+1]) for j in range(len(g))]
        for j,left,right,area,i1,i2,poly in local:
            value=geo.D(cumulative[j]);error=i2-2*value*i1+value*value*area
            if error<0:raise ArithmeticError('negative_exact_squared_error')
            if j>0:orthogonal=max(orthogonal,float(abs(i1-value*area)))
            regions.append(dict(error_squared=float(error),integral=float(value*area)))
            if len(means)<4:means.append(float(value))
    projected=geo.Field(field.domain,tuple(strips))
    info=dict(kind=kind,error_l2=math.sqrt(math.fsum(r['error_squared'] for r in regions)),
        reference_integral=math.fsum(reference),reference_l2=math.sqrt(math.fsum(reference_sq)),
        projected_integral=math.fsum(r['integral'] for r in regions),orthogonality_max=orthogonal,
        coefficients=[list(s.strengths) for s in strips])
    return projected,info

def uniform(m):
    fronts=tuple((Fraction(4*j+1,4*m),Fraction(0)) for j in range(m))
    return geo.Field(geo.rectangle(),(geo.Strip(Fraction(0),Fraction(1),fronts,tuple([0.]*m)),))

def independent_trials(native):
    records=[];quadrature=[]
    for m in [4,8,16,32]:
        h=1/m
        for kind in [1,2]:
            field,r=project(uniform(m),kind,Fraction(1));trace=math.fsum(field.strips[0].strengths)
            expected_trace=3*h/8 if kind==1 else 3*h*h/16
            record=dict(m=m,kind=kind,trace_error=1-trace,expected_trace=expected_trace,
                trace_difference=abs(1-trace-expected_trace),orthogonality=r['orthogonality_max'])
            if kind==1:record.update(error_squared_difference=abs(r['error_l2']**2-(h*h/12-11*h**3/256)),
                integral_difference=abs(r['projected_integral']-(.5-h*h/32)))
            records.append(record)
    shapes=[]
    u=uniform(4);s=u.strips[0]
    for j,left in enumerate(s.fronts):shapes.append((s,left,s.fronts[j+1] if j+1<4 else u.domain.pieces[0][3],u.domain.pieces[0]))
    for s in list(native.strips[:2])+list(native.strips[-2:]):
        piece=native.domain.piece((s.low+s.high)/2)
        for j in [0,len(s.fronts)-1]:shapes.append((s,s.fronts[j],s.fronts[j+1] if j+1<len(s.fronts) else piece[3],piece))
    monomials=[{k:Fraction(1)} for k in [(0,0),(1,0),(0,1),(2,0),(1,1),(0,2)]]
    for s,left,right,piece in shapes:
        polynomials=list(monomials)
        for kind in [1,2]:
            p=target(piece,kind,Fraction(1));area=integral(monomials[0],left,right,s.low,s.high)
            mean=integral(p,left,right,s.low,s.high)/area
            polynomials.extend([p,multiply(add(p,{(0,0):-mean}),add(p,{(0,0):-mean}))])
        differences=[abs(float(integral(p,left,right,s.low,s.high))-gauss(p,left,right,s.low,s.high)) for p in polynomials]
        quadrature.append(dict(low=float(s.low),high=float(s.high),maximum_difference=max(differences),polynomials=len(differences)))
    return dict(uniform=records,quadrature=quadrature,checks=dict(
        closed_uniform=all(max(r.get('error_squared_difference',0),r.get('integral_difference',0),r['trace_difference'])<=1e-12 for r in records),
        orthogonality=all(r['orthogonality']<=1e-12 for r in records),
        independent_moments=len(quadrature)==12 and all(r['maximum_difference']<=1e-11 for r in quadrature)))

def signed_summary(intervals):
    terms=[(float(hi-lo),d) for lo,hi,d in intervals]
    positive=math.fsum(w*max(d,0) for w,d in terms);negative=math.fsum(w*min(d,0) for w,d in terms)
    signed=math.fsum(w*d for w,d in terms);absolute=math.fsum(w*abs(d) for w,d in terms)
    return dict(positive=positive,negative=negative,signed=signed,l1=absolute,
        l2_squared=math.fsum(w*d*d for w,d in terms),cancellation_ratio=absolute/max(abs(signed),1e-12),floor_applied=abs(signed)<1e-12)

def trace_breakdown(coarse,fine):
    boundaries=[geo.D(-4+.4*j) for j in range(21)]
    z=[geo.D(x) for x in [-4,-3.6,-.4,.4,3.6,4]]
    cuts=sorted(set(boundaries+z+[y for f in [coarse,fine] for s in f.strips for y in [s.low,s.high]]))
    intervals=[];bins=[[] for _ in range(20)];zones={k:[] for k in ['left_tip','central','interior','right_tip']};rows=[]
    for lo,hi in zip(cuts[:-1],cuts[1:]):
        mid=(lo+hi)/2;d=math.fsum(fine.active(mid).strengths)-math.fsum(coarse.active(mid).strengths)
        item=(lo,hi,d);intervals.append(item)
        band=min(sum(y<=mid for y in boundaries)-1,19)
        name='left_tip' if mid<z[1] else 'right_tip' if mid>z[4] else 'central' if z[2]<mid<z[3] else 'interior'
        bins[band].append(item);zones[name].append(item)
        rows.append(dict(low=float(lo),high=float(hi),difference=d,band=band,zone=name))
    total=signed_summary(intervals);bands=[signed_summary(s) for s in bins];zone_values={k:signed_summary(s) for k,s in zones.items()}
    preservation=max(abs(math.fsum(v[k] for v in group)-total[k]) for group in [bands,list(zone_values.values())] for k in ['positive','negative','signed','l1','l2_squared'])
    return dict(total=total,bands=bands,zones=zone_values,preservation=preservation,intervals=rows)

def main():
    begin=time.monotonic();output=FOLDER/'avl-projection-trace-result.json'
    if output.exists():raise FileExistsError(output)
    parent_path=FOLDER/'avl-field-overlay-result.json';parent=json.loads(parent_path.read_text())
    inputs={};fields={};native=[];projected={};records=[];series=[];traces=[];tests=None;checks={};stop=None
    def budget():
        if time.monotonic()-begin>600:raise TimeoutError('projection_time_budget')
        if memory.memory_bytes()>512*1024**2:raise MemoryError('projection_memory_budget')
    try:
        for family,name in [('original','avl-fixed-band-grid-runs.json'),('aligned','avl-band-aligned-grid-runs.json')]:
            for run in json.loads((FOLDER/name).read_text())['runs']:
                if run['alpha']!=1:continue
                field,r=prior.load_field(run,family);key=(family,run['case'],run['m'],run['n']);fields[key]=field;native.append(r)
                for info in r['inputs'].values():inputs[info['path']]=info['sha256']
        tests=independent_trials(fields[('original','FW-V01',4,12)])
        zero=signed_summary([(Fraction(0),Fraction(1,2),1.),(Fraction(1,2),Fraction(1),-1.)])
        tests['zero_signed']=zero;tests['checks']['zero_signed_floor']=zero['signed']==0 and zero['l1']==1 and zero['floor_applied'] and zero['cancellation_ratio']==1e12
        if not all(tests['checks'].values()):raise RuntimeError('independent_projection_tests_failed')
        for key,field in fields.items():
            for kind in [1,2]:
                p,rec=project(field,kind,geo.D(prior.old.C));projected[key+(kind,)]=p
                obs,_,_=geo.overlay([p],observe=True);means=np.array(obs['cell_integrals'][0])/obs['cell_areas']
                rec.update(family=key[0],case=key[1],m=key[2],n=key[3],observations=obs,cell_means=means.tolist(),
                    relative_error_l2=rec['error_l2']/rec['reference_l2'],integral_difference=abs(obs['integrals'][0]-rec['projected_integral']))
                records.append(rec);budget()
            if len(records)%20==0:print(json.dumps(dict(stage='projection',completed=len(records),seconds=time.monotonic()-begin)),flush=True)
            write(ROOT/'build/analysis/avl-projection-trace-progress.json',dict(projections=len(records),last=key))
        lookup={(r['family'],r['case'],r['m'],r['n'],r['kind']):r for r in records}
        for s in parent['series']:
            keys=[(s['family'],s['case'],m,n) for m,n in s['meshes']]
            for j in range(3):
                breakdown=trace_breakdown(fields[keys[j]],fields[keys[j+1]])
                actual_change=s['loads']['lift']['values'][j+1]-s['loads']['lift']['values'][j]
                breakdown.update(family=s['family'],case=s['case'],direction=s['direction'],pair=j,
                    lift_change=2*breakdown['total']['signed']/prior.old.S,lift_change_difference=abs(2*breakdown['total']['signed']/prior.old.S-actual_change),
                    previous_trace_difference=max(abs(breakdown['total']['l1']-s['pairs'][j]['trace']['l1']),abs(math.sqrt(breakdown['total']['l2_squared'])-s['pairs'][j]['trace']['l2'])))
                traces.append(breakdown)
            for kind in [1,2]:
                fs=[projected[k+(kind,)] for k in keys];rr=[lookup[k+(kind,)] for k in keys]
                pairs=[];triples=[]
                for j in range(3):
                    obs,a,v=geo.overlay(fs[j:j+2]);pairs.append(dict(norms=geo.metrics(a,v[:,1]-v[:,0]),area=obs['area'],preservation=obs['integral_preservation']))
                for j in range(2):
                    _,a,v=geo.overlay(fs[j:j+3]);triples.append(prior.cosine(a,v[:,1]-v[:,0],v[:,2]-v[:,1]))
                errors=np.array([r['error_l2'] for r in rr]);orders=np.log2(errors[:-1]/errors[1:])
                cells=[np.array(r['cell_means']) for r in rr];weights=np.array(rr[-1]['observations']['cell_areas'])
                cd=[geo.metrics(weights,cells[j+1]-cells[j])['l2'] for j in range(3)]
                cc=[prior.cosine(weights,cells[j+1]-cells[j],cells[j+2]-cells[j+1]) for j in range(2)]
                series.append(dict(family=s['family'],case=s['case'],direction=s['direction'],kind=kind,
                    known_errors_l2=errors.tolist(),known_error_orders=orders.tolist(),known_errors_decrease=bool(np.all(np.diff(errors)<0)),
                    field_diagnostic=geo.diagnostic([r['norms']['l2'] for r in pairs],triples,rr[-1]['observations']['norms'][0]['l2']),
                    cell_diagnostic=geo.diagnostic(cd,cc,geo.metrics(weights,cells[-1])['l2']),pairs=pairs))
                budget()
            print(json.dumps(dict(stage='series',family=s['family'],case=s['case'],direction=s['direction'],completed=len(series),seconds=time.monotonic()-begin)),flush=True)
        checks=dict(independent_analytic_and_quadrature=all(tests['checks'].values()),eighty_projections=len(records)==80 and len(series)==24,
            projection_orthogonality=all(r['orthogonality_max']<=1e-12 for r in records),
            projected_observation_preservation=all(r['integral_difference']<=1e-10 and abs(r['observations']['area']-prior.old.S)<=1e-10 for r in records) and all(p['preservation']<=1e-10 and abs(p['area']-prior.old.S)<=1e-10 for s in series for p in s['pairs']),
            original_trace_reproduction=len(traces)==36 and all(max(r['lift_change_difference'],r['previous_trace_difference'])<=1e-11 for r in traces),
            zone_band_conservation=all(r['preservation']<=1e-10 for r in traces),inputs_unchanged=all(digest(ROOT/name)==sha for name,sha in inputs.items()),
            resource_budget=time.monotonic()-begin<600 and memory.memory_bytes()<512*1024**2)
    except Exception as exc:stop=repr(exc);raise
    finally:
        result=dict(protocol='rd/aero/avl-reference/records/avl-projection-trace-experiment.md',protocol_sha256=digest(FOLDER/'avl-projection-trace-experiment.md'),
            source_sha256=digest(Path(__file__)),dependency_sha256={Path(m.__file__).name:digest(Path(m.__file__)) for m in [geo,prior,memory]},
            parent_sha256={name:digest(FOLDER/name) for name in ['avl-field-overlay-result.json','avl-fixed-band-grid-runs.json','avl-band-aligned-grid-runs.json']},
            input_sha256=inputs,independent=tests,projections=records,projection_series=series,actual_trace=traces,
            checks=checks,methods_criteria_pass=bool(checks) and all(checks.values()),stop_reason=stop,
            reference_convergence_accepted=False,manufactured_projection_is_physical_error_estimate=False,
            seconds=time.monotonic()-begin,peak_process_working_bytes=memory.memory_bytes())
        write(output,result)
        if output.stat().st_size>64*1024**2:raise RuntimeError('projection_disk_budget')
    print(json.dumps(dict(checks=checks,seconds=result['seconds'],peak=result['peak_process_working_bytes'],bytes=output.stat().st_size),indent=2))

if __name__=='__main__':main()
