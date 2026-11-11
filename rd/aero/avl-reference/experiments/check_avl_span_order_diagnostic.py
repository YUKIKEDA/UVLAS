"""Post-hoc signed/vector diagnostics; never estimates a continuum tail or limit."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,math,time
import numpy as np
import check_avl_conservative_observation as observation
import avl_field_overlay_exact as geo
import check_q9_plate_convergence as memory
ROOT=rd_paths.ROOT;F=rd_paths.RECORDS
LEVELS=(160,320,640,1280)
MODELS={'h_h2':(1,2),'h2_h3':(2,3),'h2_h4':(2,4),'h2_h2log':(2,'log')}
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def basis(j,p):
    h=Q(1,2**j)
    return h*h*(-j) if p=='log' else h**p
def design(name):
    if name not in MODELS:raise ValueError('unknown_model')
    p,q=MODELS[name];a=[basis(j+1,p)-basis(j,p) for j in range(3)];b=[basis(j+1,q)-basis(j,q) for j in range(3)]
    det=a[0]*b[1]-a[1]*b[0]
    if det==0:raise ValueError('singular_model')
    inverse=((b[1]/det,-b[0]/det),(-a[1]/det,a[0]/det))
    weights=(a[2]*inverse[0][0]+b[2]*inverse[1][0],a[2]*inverse[0][1]+b[2]*inverse[1][1])
    return a,b,inverse,weights
def checked(values,levels=LEVELS):
    if tuple(levels)!=LEVELS:raise ValueError('unsupported_levels')
    v=np.asarray(values,dtype=float)
    if v.ndim==1:v=v[:,None]
    if v.ndim!=2 or v.shape[0]!=4 or v.shape[1]==0:raise ValueError('shape')
    if not np.all(np.isfinite(v)):raise ValueError('nonfinite')
    return v
def cosine(a,b):
    den=float(np.linalg.norm(a)*np.linalg.norm(b));return float(a@b/den) if den else None
def ratio(value,scale):return float(value/scale) if scale>0 else None
def signed(values):
    pos=math.fsum(float(v) for v in values if v>0);neg=math.fsum(float(v) for v in values if v<0);total=math.fsum(float(v) for v in values);absolute=pos-neg
    return dict(total=total,positive=pos,negative=neg,absolute_sum=absolute,cancellation_ratio=abs(total)/absolute if absolute else None)
def fit(values,name,levels=LEVELS):
    v=checked(values,levels);d=np.diff(v,axis=0);a,b,inv,w=design(name)
    ca=float(inv[0][0])*d[0]+float(inv[0][1])*d[1];cb=float(inv[1][0])*d[0]+float(inv[1][1])*d[1]
    terms=np.array([float(a[2])*ca,float(b[2])*cb]);pred=float(w[0])*d[0]+float(w[1])*d[1]
    matrix=np.array([[float(a[j]),float(b[j])] for j in range(2)])
    alternative=np.array([float(a[2]),float(b[2])])@np.linalg.solve(matrix,d[:2])
    norm=float(np.linalg.norm(d[2]));termscale=sum(float(np.linalg.norm(t)) for t in terms)
    return dict(model=name,rational_weights=[str(x) for x in w],design=matrix.tolist(),design_condition=float(np.linalg.cond(matrix)),
        predicted_last_difference=pred.tolist(),actual_last_difference=d[2].tolist(),residual=(pred-d[2]).tolist(),
        relative_residual=ratio(np.linalg.norm(pred-d[2]),norm),prediction_cosine=cosine(pred,d[2]),
        last_terms=terms.tolist(),model_cancellation_ratio=ratio(np.linalg.norm(pred),termscale),term_cosine=cosine(*terms),
        independent_prediction_difference=float(np.max(abs(pred-alternative))),
        coefficient_a=ca.tolist(),coefficient_b=cb.tolist(),scope='post-hoc model diagnostic; no limit or tail bound')
def baseline(values):
    d=np.diff(checked(values),axis=0);out=[]
    for p in [1.,2.,3.,4.]:
        pred=d[0]*2**(-2*p);out.append(dict(order=p,relative_last_residual=ratio(np.linalg.norm(pred-d[2]),np.linalg.norm(d[2]))))
    n0=np.linalg.norm(d[0]);n1=np.linalg.norm(d[1]);p=math.log2(n0/n1) if min(n0,n1)>0 else None
    pred=d[1]*(n1/n0) if p is not None else None
    return dict(fixed=out,first_pair_order=p,constant_order_last_residual=ratio(np.linalg.norm(pred-d[2]),np.linalg.norm(d[2])) if pred is not None else None)
def manufacturing():
    records=[];a=np.array([1.,-2.,.5]);b=np.array([-.3,.7,1.2]);offset=np.array([3.,4.,-1.])
    for name,(p,q) in MODELS.items():
        values=np.array([offset+float(basis(j,p))*a+float(basis(j,q))*b for j in range(4)])
        r=fit(values,name);error=float(np.linalg.norm(r['residual']));assert error<1e-12
        records.append(dict(kind=name,prediction_error=error,passed=True))
    # Exact analytic example: the first three differences have the same direction,
    # apparent order grows like the data, and the next difference reverses sign.
    values=[Q(5,2)*Q(1,2**j)-Q(26,3)*Q(1,2**j)**2 for j in range(5)]
    dif=[values[j+1]-values[j] for j in range(4)];assert dif==[Q(21,4),Q(1),Q(3,32),Q(-7,128)]
    records.append(dict(kind='competing_terms_sign_reversal',exact_differences=[str(x) for x in dif],
        orders=[math.log2(float(abs(dif[j]/dif[j+1]))) for j in range(3)],passed=True))
    three=np.array([offset+float(basis(j,1))*a+float(basis(j,2))*b+float(basis(j,3))*np.array([.4,.1,-.2]) for j in range(4)])
    error=np.linalg.norm(fit(three,'h_h2')['residual']);assert error>1e-4
    records.append(dict(kind='unmodeled_third_term',prediction_error=float(error),passed=True))
    c=signed([1.,-1.]);assert c['total']==0 and np.linalg.norm([1.,-1.])>1
    records.append(dict(kind='scalar_hides_distribution',signed=c,norm=math.sqrt(2),passed=True))
    constant=np.tile(offset,(4,1));r=fit(constant,'h_h2');assert r['relative_residual'] is None
    records.append(dict(kind='constant_no_relative_claim',passed=True))
    roots=[Q(1,2**j) for j in range(4)]
    def ambiguity(h):return 64*math.prod(h-x for x in roots)
    assert all(ambiguity(x)==0 for x in roots) and ambiguity(Q(0))==1
    records.append(dict(kind='finite_samples_do_not_determine_limit',sample_values=[str(ambiguity(x)) for x in roots],limit_difference=1,passed=True))
    failures=[]
    for name,values,model,levels in [('nan',[[0],[1],[2],[float('nan')]],'h_h2',LEVELS),('shape',[[0],[1],[2]],'h_h2',LEVELS),
        ('levels',np.zeros((4,1)),'h_h2',(160,320,640,640)),('model',np.zeros((4,1)),'other',LEVELS)]:
        before=repr(values)
        try:fit(values,model,levels);refused=False
        except ValueError:refused=True
        assert refused and repr(values)==before;failures.append(dict(kind=name,refused=True,input_unchanged=True))
    return dict(cases=records,refusals=failures,all_criteria_pass=True)
def regions(bands):
    centers=np.array([b['center'] for b in bands]);return {'root':abs(centers)<=.8,'tip':abs(centers)>=3.2,'interior':(abs(centers)>.8)&(abs(centers)<3.2)}
def moment_parts(bands,case):
    domain=geo.wing_domain(case);parts=[]
    for b in bands:
        y=geo.D(b['center']);te=domain.piece(y)[3];weight=float(geo.at(te,y)-observation.OX);slope=float(te[1])
        parts.append([-2*(weight*b['T0']+slope*b['T1'])/(observation.S*observation.C),2*b['K0']/(observation.S*observation.C),b['Cy']/(observation.S*observation.C)])
    return np.array(parts)
def main():
    output=F/'avl-span-order-diagnostic-result.json'
    if output.exists():raise FileExistsError(output)
    begin=time.monotonic();parent_path=F/'avl-span1280-observation-result.json';parent=json.loads(parent_path.read_text());assert parent['methods_criteria_pass']
    unit=manufacturing();channels=[];scalars=[];checks={};maximum_reproduction=0.;maximum_moment=0.;maximum_partition=0.;maximum_fit=0.
    for case in ['FW-V01','FW-V02']:
        meshes=[next(m for m in parent['meshes'] if (m['case'],m['m'],m['n'])==(case,24,n)) for n in LEVELS]
        bands=[m['bands'] for m in meshes];vectors=[observation.channels(b) for b in bands];masks=regions(bands[0]);counts={k:int(np.sum(v)) for k,v in masks.items()};assert sum(counts.values())==20
        saved=next(s for s in parent['series'] if s['case']==case and s['label']=='span1280')
        for name in vectors[0]:
            v=np.array([x[name] for x in vectors]);d=np.diff(v,axis=0);diagnostic=observation.norm_diagnostic(list(v));old=saved['channels'][name]
            maximum_reproduction=max(maximum_reproduction,max(abs(np.array(diagnostic['differences'])-np.array(old['differences']))),max(abs(np.array(diagnostic['orders'])-np.array(old['orders']))))
            stages=[]
            for j in range(3):
                energy=d[j,:20]**2+(d[j,20:]**2 if name=='tau_p1' else 0);norm2=float(d[j]@d[j]);regional={r:float(np.sum(energy[mask])) for r,mask in masks.items()}
                maximum_partition=max(maximum_partition,abs(sum(regional.values())-norm2))
                stages.append(dict(pair=[LEVELS[j],LEVELS[j+1]],difference=d[j].tolist(),norm=math.sqrt(norm2),band_energy=energy.tolist(),regional_energy=regional,
                    regional_fractions={r:ratio(e,norm2) for r,e in regional.items()}))
            projections=[]
            for j in range(2):
                factor=float(d[j+1]@d[j]/(d[j]@d[j]));perp=d[j+1]-factor*d[j]
                projections.append(dict(pair=j,factor=factor,orthogonal_fraction=ratio(np.linalg.norm(perp),np.linalg.norm(d[j+1])),cosine=cosine(d[j],d[j+1])))
            local=[];floor=1e-12*diagnostic['scale']
            for index in range(v.shape[1]):
                ds=d[:,index];orders=[]
                for j in range(2):orders.append(dict(state='small' if min(abs(ds[j]),abs(ds[j+1]))<=floor else ('sign_reversal' if ds[j]*ds[j+1]<0 else 'same_sign'),
                    order=math.log2(abs(ds[j]/ds[j+1])) if min(abs(ds[j]),abs(ds[j+1]))>floor and ds[j]*ds[j+1]>0 else None))
                local.append(dict(component=index,band=index%20,part='p1' if index>=20 else 'mean',differences=ds.tolist(),orders=orders))
            fits=[fit(v,m) for m in MODELS];maximum_fit=max(maximum_fit,max(r['independent_prediction_difference'] for r in fits))
            channels.append(dict(case=case,channel=name,values=v.tolist(),diagnostic=diagnostic,region_counts=counts,stages=stages,projections=projections,local=local,models=fits,baselines=baseline(v)))
        parts=np.array([moment_parts(b,case) for b in bands]);direct=np.array([[b['moment_point'] for b in bs] for bs in bands]);maximum_moment=max(maximum_moment,float(np.max(abs(np.sum(parts,axis=2)-direct))))
        for name in ['lift','moment_line','moment_point']:
            values=np.array([[b[name] for b in bs] for bs in bands]);totals=np.array([math.fsum(row) for row in values]);diff=np.diff(values,axis=0);records=[]
            for j in range(3):
                s=signed(diff[j]);regional={r:signed(diff[j,mask]) for r,mask in masks.items()};maximum_partition=max(maximum_partition,abs(sum(x['total'] for x in regional.values())-s['total']))
                item=dict(pair=[LEVELS[j],LEVELS[j+1]],bands=diff[j].tolist(),signed=s,regions=regional)
                if name=='moment_point':
                    pd=parts[j+1]-parts[j];item['physical_parts']={label:signed(pd[:,i]) for i,label in enumerate(['tau_position','k','couple'])}
                    maximum_moment=max(maximum_moment,abs(sum(x['total'] for x in item['physical_parts'].values())-s['total']))
                records.append(item)
            scalars.append(dict(case=case,observable=name,values=totals.tolist(),stages=records,models=[fit(totals,m) for m in MODELS],baselines=baseline(totals)))
    checks=dict(manufacturing=unit['all_criteria_pass'],ten_channels=len(channels)==10,six_scalar_series=len(scalars)==6,
        original_diagnostic_reproduction=maximum_reproduction<=1e-12,physical_moment_reconstruction=maximum_moment<=1e-12,
        spatial_partition=maximum_partition<=1e-12,independent_model_prediction=maximum_fit<=1e-12,
        resource_budget=time.monotonic()-begin<120 and memory.memory_bytes()<512*1024**2)
    result=dict(protocol='rd/aero/avl-reference/records/avl-span-order-diagnostic-experiment.md',protocol_sha256=digest(F/'avl-span-order-diagnostic-experiment.md'),
        source_sha256=digest(__file__),dependency_sha256={Path(m.__file__).name:digest(m.__file__) for m in [observation,geo,memory]},
        parent_sha256=digest(parent_path),manufacturing=unit,channels=channels,scalars=scalars,checks=checks,
        maximum_reproduction_difference=maximum_reproduction,maximum_moment_difference=maximum_moment,maximum_partition_difference=maximum_partition,
        maximum_independent_fit_difference=maximum_fit,seconds=time.monotonic()-begin,peak_bytes=memory.memory_bytes(),
        all_criteria_pass=all(checks.values()),new_avl_runs=0,continuum_tail_estimate=None,reference_acceptance=False)
    observation.write(output,result);assert all(checks.values()),checks
    print(json.dumps(dict(checks=checks,manufactured_cases=len(unit['cases']),refusals=len(unit['refusals']),seconds=result['seconds'],
        vector_models=[dict(case=r['case'],channel=r['channel'],models={m['model']:m['relative_residual'] for m in r['models']}) for r in channels]),indent=2))
if __name__=='__main__':main()
