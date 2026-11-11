"""Predeclared reference-only tau/K moment-preserving observation experiment."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
from fractions import Fraction
from bisect import bisect_right
import hashlib,json,math,time
import numpy as np
import avl_field_overlay_exact as geo
import check_avl_field_overlay as prior
import check_avl_linear_load_wrench as wrench
import check_q9_plate_convergence as memory

ROOT=rd_paths.ROOT;FOLDER=rd_paths.RECORDS
S=prior.old.S;C=prior.old.C;OX=geo.D(.375)
BOUNDS=tuple(geo.D(-4+.4*j) for j in range(21))
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,x):p.write_text(json.dumps(x,indent=2,default=prior.scalar,allow_nan=False)+'\n')
def ip(a,lo,hi):return a[0]*(hi-lo)+a[1]*(hi*hi-lo*lo)/2
def ip2(a,lo,hi):return a[0]*a[0]*(hi-lo)+a[0]*a[1]*(hi*hi-lo*lo)+a[1]*a[1]*(hi**3-lo**3)/3

def profile(field,couples):
    if field.basis!='horseshoe':raise ValueError('basis_mismatch')
    if len(couples)!=len(field.strips):raise ValueError('couple_shape')
    if not all(math.isfinite(float(c)) for c in couples):raise ValueError('nonfinite_couple')
    rows=[]
    for s,cy in zip(field.strips,couples):
        if any(s.low<p[0]<s.high for p in field.domain.pieces):raise ValueError('unresolved_source_break')
        te=field.domain.piece((s.low+s.high)/2)[3];g=[geo.D(v) for v in s.strengths]
        tau=sum(g,Fraction(0))
        k=tuple(sum((a*(te[j]-l[j]) for a,l in zip(g,s.fronts)),Fraction(0)) for j in range(2))
        rows.append((s.low,s.high,tau,k,geo.D(cy)))
    return dict(field=field,rows=tuple(rows))

def active(model,y):return model['rows'][bisect_right([r[1] for r in model['rows']],y)]

def observe(model,bounds=BOUNDS):
    field=model['field'];bounds=tuple(geo.D(v) for v in bounds)
    if len(bounds)<2 or bounds[0]!=field.domain.low or bounds[-1]!=field.domain.high or any(b<=a for a,b in zip(bounds[:-1],bounds[1:])):raise ValueError('invalid_observation_bounds')
    bands=[]
    for lo,hi in zip(bounds[:-1],bounds[1:]):
        if any(lo<p[0]<hi for p in field.domain.pieces):raise ValueError('unresolved_weight_breakpoint')
        yc=(lo+hi)/2;h=hi-lo;te=field.domain.piece(yc)[3]
        cuts=sorted({lo,hi}|{y for r in model['rows'] for y in r[:2] if lo<y<hi})
        t0=Fraction(0);t1=Fraction(0);k0=Fraction(0);cy=Fraction(0)
        for a,b in zip(cuts[:-1],cuts[1:]):
            _,_,tau,k,c=active(model,(a+b)/2)
            t0+=tau*(b-a);t1+=tau*((b-yc)**2-(a-yc)**2)/2
            k0+=ip(k,a,b);cy+=c*(b-a)
        j=(geo.at(te,yc)-OX)*t0+te[1]*t1-k0
        bands.append(dict(low=lo,high=hi,width=h,center=yc,T0=t0,T1=t1,K0=k0,Cy=cy,
            tau0=t0/h,tau1=6*t1/h**2,kbar=k0/h,jbar=j/h,jcbar=(j-cy/2)/h,J=j,
            lift=2*float(t0)/S,moment_line=-2*float(j)/(S*C),moment_point=float(-2*j+cy)/(S*C)))
    return bands

def public(bands):return [{k:float(v) for k,v in b.items()} for b in bands]
def totals(bands):return {k:math.fsum(b[k] for b in bands) for k in ['lift','moment_line','moment_point']}
def channels(bands):
    h=np.array([float(b['width']) for b in bands]);root=np.sqrt(h)
    return dict(tau0=root*np.array([float(b['tau0']) for b in bands]),
        tau_p1=np.concatenate([root*np.array([float(b['tau0']) for b in bands]),np.sqrt(h/3)*np.array([float(b['tau1']) for b in bands])]),
        k0=root*np.array([float(b['kbar']) for b in bands]),j_line=root*np.array([float(b['jbar']) for b in bands]),j_point=root*np.array([float(b['jcbar']) for b in bands]))

def raw_norms(a,b):
    if a['field'].domain!=b['field'].domain:raise ValueError('domain_mismatch')
    domain=a['field'].domain
    ys=sorted({y for m in [a,b] for r in m['rows'] for y in r[:2]}|{p[0] for p in domain.pieces})
    tau=[];k=[]
    for lo,hi in zip(ys[:-1],ys[1:]):
        ra=active(a,(lo+hi)/2);rb=active(b,(lo+hi)/2)
        dt=rb[2]-ra[2];dk=tuple(y-x for x,y in zip(ra[3],rb[3]))
        tau.append(float(dt*dt*(hi-lo)));k.append(float(ip2(dk,lo,hi)))
    w=math.fsum(float(ip2((p[3][0]-OX,p[3][1]),p[0],p[1])) for p in domain.pieces)
    return math.sqrt(math.fsum(tau)),math.sqrt(math.fsum(k)),math.sqrt(w),float(domain.high-domain.low)

def pair_bounds(a,b,oa,ob):
    dt,dk,w,breadth=raw_norms(a,b);va=channels(oa);vb=channels(ob);ta=totals(oa);tb=totals(ob)
    norms={k:float(np.linalg.norm(vb[k]-va[k])) for k in va}
    direct={k:abs(tb[k]-ta[k]) for k in ta}
    delta_couple=abs((tb['moment_point']-tb['moment_line'])-(ta['moment_point']-ta['moment_line']))
    line=2/(S*C)*(w*dt+math.sqrt(breadth)*dk)
    projected_line=2/(S*C)*(w*norms['tau_p1']+math.sqrt(breadth)*norms['k0'])
    bounds=dict(raw_lift=2*math.sqrt(breadth)/S*dt,raw_line=line,raw_point=line+delta_couple,
        projected_lift=2*math.sqrt(breadth)/S*norms['tau_p1'],projected_line=projected_line,
        projected_point=projected_line+delta_couple,
        moment_density_line=2*math.sqrt(breadth)/(S*C)*norms['j_line'],
        moment_density_point=2*math.sqrt(breadth)/(S*C)*norms['j_point'])
    margins={k:bounds[k]-direct[name] for k,name in [('raw_lift','lift'),('projected_lift','lift'),
        ('raw_line','moment_line'),('projected_line','moment_line'),('moment_density_line','moment_line'),
        ('raw_point','moment_point'),('projected_point','moment_point'),('moment_density_point','moment_point')]}
    return dict(actual_absolute_changes=direct,bounds=bounds,margins=margins,
        raw_norms=dict(tau=dt,k=dk,weight=w),projected_norms=norms,
        contraction_excess=max(norms['tau_p1']**2-dt**2,norms['k0']**2-dk**2,0.),
        bound_scope='difference of these two stored fields; not continuum error')

def point_couples(field,theta):
    densities=[];point_terms=[]
    for s in field.strips:
        y=s.low+(s.high-s.low)*geo.D(theta);mid=(s.low+s.high)/2
        density=sum((-2*geo.D(g)*(geo.at(l,y)-geo.at(l,mid)) for l,g in zip(s.fronts,s.strengths)),Fraction(0))
        densities.append(density)
        point_terms.append(float(-2*sum((geo.D(g)*(geo.at(l,y)-OX) for l,g in zip(s.fronts,s.strengths)),Fraction(0))*(s.high-s.low))/(S*C))
    return tuple(densities),math.fsum(point_terms)

def independent_gauss(field,theta,bounds):
    nodes,weights=np.polynomial.legendre.leggauss(4);result=[]
    for lo,hi in zip(bounds[:-1],bounds[1:]):
        yc=float((lo+hi)/2);cuts=sorted({lo,hi}|{y for s in field.strips for y in [s.low,s.high] if lo<y<hi})
        terms=[]
        for low,high in zip(cuts[:-1],cuts[1:]):
            s=field.active((low+high)/2);piece=field.domain.piece((low+high)/2)
            # Separate point-force/couple computation in float64 vector algebra.
            couple=[]
            for l,g in zip(s.fronts,s.strengths):
                aa=np.array([float(geo.at(l,s.low)),float(s.low),0.]);bb=np.array([float(geo.at(l,s.high)),float(s.high),0.])
                p=aa+theta*(bb-aa);force=2*g*np.cross(np.array([1.,0.,0.]),bb-aa)
                couple.append(float(np.cross(p-(aa+bb)/2,force)[1]))
            cy=math.fsum(couple)/float(s.high-s.low)
            for q,w in zip(nodes,weights):
                y=float(low)+(q+1)*float(high-low)/2
                le=float(piece[2][0])+float(piece[2][1])*y;te=float(piece[3][0])+float(piece[3][1])*y
                fronts=[float(l[0])+float(l[1])*y for l in s.fronts];xs=[le]+fronts+[te]
                kval=[]
                for xl,xr in zip(xs[:-1],xs[1:]):
                    samples=xl+(nodes+1)*(xr-xl)/2
                    values=np.array([math.fsum(g for g,x0 in zip(s.strengths,fronts) if x>=x0) for x in samples])
                    kval.append((xr-xl)/2*float(weights@values))
                tau=math.fsum(s.strengths);k=math.fsum(kval);factor=float(w)*float(high-low)/2
                terms.append(factor*np.array([tau,(y-yc)*tau,k,cy,(te-.375)*tau-k]))
        result.append([math.fsum(v[j] for v in terms) for j in range(5)])
    return np.array(result)

def manufacturing():
    D=geo.D;domain=geo.Domain('moment-manufactured',((Fraction(0),Fraction(1),(Fraction(0),Fraction(1,10)),(Fraction(1),Fraction(2,5))),))
    limits=[Fraction(0),Fraction(1,3),Fraction(3,5),Fraction(1)]
    fronts=((Fraction(1,5),Fraction(4,25)),(Fraction(3,5),Fraction(7,25)))
    strengths=[(1.,.5),(.5,.25),(.75,.75)];bounds=(Fraction(0),Fraction(1,5),Fraction(11,20),Fraction(9,10),Fraction(1))
    cases=[];models=[];observed=[]
    for name,factor,theta in [('zero_couple',1.,.5),('quarter_point',1.,.25),('negative',-.7,.25),('scaled',2.3,.25)]:
        field=geo.Field(domain,tuple(geo.Strip(lo,hi,fronts,tuple(factor*g for g in gs)) for lo,hi,gs in zip(limits[:-1],limits[1:],strengths)))
        cp,direct=point_couples(field,theta);model=profile(field,cp);obs=observe(model,bounds);quad=independent_gauss(field,theta,bounds)
        exact=np.array([[float(o[k]) for k in ['T0','T1','K0','Cy','J']] for o in obs])
        total=totals(obs);drop_t1=math.fsum(-2*float((geo.at(domain.pieces[0][3],o['center'])-OX)*o['T0']-o['K0'])/(S*C)+float(o['Cy'])/(S*C) for o in obs)
        cases.append(dict(name=name,maximum_gauss_difference=float(np.max(abs(exact-quad))),point_difference=abs(total['moment_point']-direct),
            omit_T1_difference=abs(total['moment_point']-drop_t1),omit_couple_difference=abs(total['moment_point']-total['moment_line']),bands=public(obs)))
        models.append(model);observed.append(obs)
    pairs=[pair_bounds(models[0],m,observed[0],o) for m,o in zip(models[1:],observed[1:])]
    field=models[0]['field'];before=hashlib.sha256(repr(field).encode()).hexdigest();failures=[]
    wing=geo.wing_domain('fold-test');half=geo.Field(wing,(geo.Strip(D(-4),D(0),((D(.75),D(-.6)),),(1.,)),geo.Strip(D(0),D(4),((D(.75),D(.6)),),(1.,))))
    actions=[('fold_weight',lambda:observe(profile(half,(0.,0.)),(D(-4),D(4)))),
        ('degenerate_band',lambda:observe(models[0],(D(0),D(.5),D(.5),D(1)))),
        ('reversed_band',lambda:observe(models[0],(D(0),D(.8),D(.2),D(1)))),
        ('nonfinite_couple',lambda:profile(field,(0.,float('nan'),0.))),
        ('couple_shape',lambda:profile(field,(0.,))),
        ('basis_mismatch',lambda:profile(geo.Field(field.domain,field.strips,'ring'),(0.,0.,0.))),
        ('domain_mismatch',lambda:raw_norms(models[0],profile(half,(0.,0.))))]
    for name,fn in actions:
        try:fn();refused=False
        except ValueError:refused=True
        failures.append(dict(kind=name,refused=refused,input_unchanged=hashlib.sha256(repr(field).encode()).hexdigest()==before))
    return dict(cases=cases,pairs=pairs,failures=failures,checks=dict(
        independent_gauss=all(c['maximum_gauss_difference']<=1e-11 for c in cases),direct_point=all(c['point_difference']<=1e-11 for c in cases),
        negative_T1=max(c['omit_T1_difference'] for c in cases)>1e-5,negative_couple=max(c['omit_couple_difference'] for c in cases)>1e-5,
        stored_field_bounds=all(min(p['margins'].values())>=-1e-10 and p['contraction_excess']<=1e-10 for p in pairs),
        failures=all(f['refused'] and f['input_unchanged'] for f in failures)))

def native_model(run,family):
    field,record=prior.load_field(run,family);folder=ROOT/run['directory']
    rows=np.loadtxt(folder/'uv-rows.txt',skiprows=2);forces=np.loadtxt(folder/'uv-forces.txt',skiprows=2);ff=forces[forces[:,2]==0]
    assert np.array_equal(ff[:,0],np.arange(1,len(rows)+1))
    a=rows[:,7:10];b=rows[:,10:13];p=ff[:,3:6];g=rows[:,13]/math.sin(math.pi/180)
    f=2*g[:,None]*np.cross(np.array([1.,0.,0.]),b-a);cp=wrench.convert(a,b,p,f)
    densities=[]
    for strip,s in zip(sorted(set(ff[:,1])),field.strips):
        values=cp[ff[:,1]==strip,1];densities.append(sum((geo.D(x) for x in values),Fraction(0))/(s.high-s.low))
    model=profile(field,tuple(densities));observed=observe(model);reference=wrench.observe(a,b,f,cp)
    total=totals(observed)
    record.update(bands=public(observed),observation_id='tau-k-band-moments-v1',
        load_identity_difference=max(abs(total['lift']-record['lift']),abs(total['moment_line']-record['moment_line']),abs(total['moment_point']-record['moment_point'])),
        independent_band_difference=max(float(np.max(abs(np.array([o['lift'] for o in observed])-reference['force'][:,2]/S))),
            float(np.max(abs(np.array([o['moment_point'] for o in observed])-reference['moment'][:,1]/(S*C))))),
        p1_moment_preservation=max(max(abs(float(o['tau0']*o['width']-o['T0'])),abs(float(o['tau1']*o['width']**2/6-o['T1']))) for o in observed))
    return model,observed,record

def norm_diagnostic(values):
    differences=[values[j+1]-values[j] for j in range(3)];magnitudes=[float(np.linalg.norm(v)) for v in differences]
    cosines=[None if magnitudes[j]*magnitudes[j+1]==0 else float(differences[j]@differences[j+1]/(magnitudes[j]*magnitudes[j+1])) for j in range(2)]
    scale=float(np.linalg.norm(values[-1]));d=geo.diagnostic(magnitudes,cosines,scale)
    d['relative_last_difference']=magnitudes[-1]/scale if scale>0 else None
    d['conditional_extrapolation']=1.25*magnitudes[-1]/(2**d['orders'][-1]-1) if d['status']=='asymptotic_candidate' else (4*max(magnitudes) if d['status']=='roundoff_candidate' else None)
    d['estimate_kind']='conditional extrapolation; not certified continuum bound'
    d['scale']=scale;return d

def main():
    begin=time.monotonic();output=FOLDER/'avl-conservative-observation-result.json'
    if output.exists():raise FileExistsError(output)
    parent=json.loads((FOLDER/'avl-field-overlay-result.json').read_text());unit=None;meshes=[];pairs=[];series=[];budgets=[];checks={};stop=None;models={};obs={}
    def budget():
        if time.monotonic()-begin>300:raise TimeoutError('observation_time_budget')
        if memory.memory_bytes()>512*1024**2:raise MemoryError('observation_memory_budget')
    try:
        unit=manufacturing()
        if not all(unit['checks'].values()):raise RuntimeError('manufactured_observation_failed')
        for family,name in [('original','avl-fixed-band-grid-runs.json'),('aligned','avl-band-aligned-grid-runs.json')]:
            for run in json.loads((FOLDER/name).read_text())['runs']:
                if run['alpha']!=1:continue
                key=(family,run['case'],run['m'],run['n']);models[key],obs[key],record=native_model(run,family)
                previous=next(r for r in parent['meshes'] if all(r[k]==record[k] for k in ['family','case','m','n']))
                record['previous_ledger_difference']=max(abs(totals(obs[key])[a]-previous[b]) for a,b in [('lift','lift'),('moment_line','moment_line'),('moment_point','moment_point')])
                meshes.append(record);budget()
            print(json.dumps(dict(stage='meshes',family=family,completed=len(meshes),seconds=time.monotonic()-begin)),flush=True)
        for s in parent['series']:
            keys=[(s['family'],s['case'],m,n) for m,n in s['meshes']]
            vectors=[channels(obs[k]) for k in keys]
            series.append(dict(family=s['family'],case=s['case'],direction=s['direction'],channels={k:norm_diagnostic([v[k] for v in vectors]) for k in vectors[0]}))
            for j in range(3):
                p=pair_bounds(models[keys[j]],models[keys[j+1]],obs[keys[j]],obs[keys[j+1]])
                p.update(family=s['family'],case=s['case'],direction=s['direction'],pair=j);pairs.append(p)
            budget()
        for family in ['original','aligned']:
            for case in ['FW-V01','FW-V02']:
                channels_budget={}
                for channel in ['tau0','tau_p1','k0','j_line','j_point']:
                    axes={d:next(s for s in series if (s['family'],s['case'],s['direction'])==(family,case,d))['channels'][channel] for d in ['chord','span']}
                    estimates=[v['conditional_extrapolation'] for v in axes.values()]
                    available=all(v is not None for v in estimates);scale=axes['chord']['scale']
                    channels_budget[channel]=dict(axes={k:v['status'] for k,v in axes.items()},
                        combined_conditional_relative=sum(estimates)/scale if available and scale>0 else None,
                        state='conditional_estimate_only' if available else 'unassessed_missing_axis',
                        missing_axes=[k for k,v in axes.items() if v['conditional_extrapolation'] is None])
                budgets.append(dict(family=family,case=case,channels=channels_budget,reference_acceptance_eligible=False))
        checks=dict(manufactured=all(unit['checks'].values()),forty_meshes=len(meshes)==40 and len(series)==12,
            load_identities=all(max(r['load_identity_difference'],r['previous_ledger_difference'])<=1e-10 for r in meshes),
            independent_bands=all(r['independent_band_difference']<=1e-10 for r in meshes),
            p1_moment_preservation=all(r['p1_moment_preservation']<=1e-10 for r in meshes),
            stored_field_upper_bounds=len(pairs)==36 and all(min(r['margins'].values())>=-1e-10 for r in pairs),
            projection_contraction=all(r['contraction_excess']<=1e-10 for r in pairs),
            inputs_unchanged=all(digest(ROOT/v['path'])==v['sha256'] for r in meshes for v in r['inputs'].values()),
            resource_budget=time.monotonic()-begin<300 and memory.memory_bytes()<512*1024**2)
    except Exception as exc:stop=repr(exc);raise
    finally:
        result=dict(protocol='rd/aero/avl-reference/records/avl-conservative-observation-experiment.md',protocol_sha256=digest(FOLDER/'avl-conservative-observation-experiment.md'),
            source_sha256=digest(Path(__file__)),dependency_sha256={Path(m.__file__).name:digest(Path(m.__file__)) for m in [geo,prior,wrench,memory]},
            parent_sha256={name:digest(FOLDER/name) for name in ['avl-field-overlay-result.json','avl-fixed-band-grid-runs.json','avl-band-aligned-grid-runs.json']},
            manufactured=unit,meshes=meshes,pairs=pairs,series=series,budgets=budgets,checks=checks,
            methods_criteria_pass=bool(checks) and all(checks.values()),stop_reason=stop,reference_convergence_accepted=False,
            continuum_error_bound_available=False,seconds=time.monotonic()-begin,peak_process_working_bytes=memory.memory_bytes())
        write(output,result)
        if output.stat().st_size>64*1024**2:raise RuntimeError('observation_disk_budget')
    print(json.dumps(dict(checks=checks,seconds=result['seconds'],peak=result['peak_process_working_bytes'],bytes=output.stat().st_size),indent=2))

if __name__=='__main__':main()
